#include <mpi.h>
#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

// RSA MPI
#include <RSAMPI/fields.h>

#include <RSAMPI/exanb_naive_resolve.hxx>
#include <radius_generator.hxx>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbResolveCandidatesPass : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"Main grid - already-committed spheres"});
  ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{"Thickness (in cells) of the grid's ghost shell"});
  ADD_SLOT(GridT, candidate_grid, INPUT_OUTPUT,
           DocString{"Grid holding every rank's still-alive candidates (real=mine, ghost=neighbors')"});
  ADD_SLOT(bool, any_changed, INPUT_OUTPUT, true, DocString{"True as long as some rank eliminated a candidate"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius generator state"});
  ADD_SLOT(bool, log_undecided, INPUT, false,
           DocString{"Log the (total / avg / max per rank) count of not-yet-confirmed candidates each pass"});

  inline void execute() override final {
    const auto stats =
        exanb_naive::resolve_candidates_pass(*grid, *candidate_grid, *ghost_layer, RSARadiusGenerator->get_max_radius());
    candidate_grid->rebuild_particle_offsets();

    int local_changed = stats.changed ? 1 : 0;
    int global_changed = 0;
    MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    int64_t total_undecided = 0;
    MPI_Allreduce(&stats.undecided, &total_undecided, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    // no undecided candidate left anywhere -> nothing left a future pass could
    // do, even if this pass itself changed something.
    *any_changed = (global_changed != 0) && (total_undecided != 0);

    if (*log_undecided) {
      int size = 1, rank = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      int64_t max_undecided = 0;
      MPI_Allreduce(&stats.undecided, &max_undecided, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
      if (rank == 0) {
        onika::lout << "Undecided candidates: total=" << total_undecided << " avg="
                    << (double(total_undecided) / size) << " max=" << max_undecided << std::endl;
      }
    }
  }
};

template <class GridT>
using RSAMPIExanbResolveCandidatesPassTmpl = RSAMPIExanbResolveCandidatesPass<GridT>;

ONIKA_AUTORUN_INIT(exanb_resolve_candidates_pass) {
  OperatorNodeFactory::instance()->register_factory(
      "exanb_resolve_candidates_pass", exanb::make_grid_variant_operator<RSAMPIExanbResolveCandidatesPassTmpl>);
}
}  // namespace rsa_mpi
