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

#include <RSAMPI/exanb_naive_algorithm.hxx>

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

  inline void execute() override final {
    const bool local_changed = exanb_naive::resolve_candidates_pass(*grid, *candidate_grid, *ghost_layer);
    candidate_grid->rebuild_particle_offsets();

    int local = local_changed ? 1 : 0;
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    *any_changed = (global != 0);
  }
};

template <class GridT>
using RSAMPIExanbResolveCandidatesPassTmpl = RSAMPIExanbResolveCandidatesPass<GridT>;

ONIKA_AUTORUN_INIT(exanb_resolve_candidates_pass) {
  OperatorNodeFactory::instance()->register_factory(
      "exanb_resolve_candidates_pass", exanb::make_grid_variant_operator<RSAMPIExanbResolveCandidatesPassTmpl>);
}
}  // namespace rsa_mpi
