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

#include <RSAMPI/exanb_naive_grid.hxx>
#include <RSAMPI/exanb_naive_types.hxx>
#include <limits>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbCandidatesToGrid : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"exaNBody grid (used only for its geometry)"});
  ADD_SLOT(exanb_naive::Candidates, candidates, INPUT, REQUIRED, DocString{"This round's local candidates"});
  ADD_SLOT(GridT, candidate_grid, INPUT_OUTPUT,
           DocString{"Freshly (re)built grid holding this round's local candidates"});
  ADD_SLOT(bool, any_changed, INPUT_OUTPUT, true, DocString{"Reset each round for resolve_candidates_step"});

  inline void execute() override final {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    *candidate_grid = exanb_naive::make_scratch_grid(*grid);
    *any_changed = true;

    const double cs = grid->cell_size();
    const exanb::Vec3d origin = grid->origin();
    const exanb::IJK off = grid->offset();
    const ssize_t gl = static_cast<ssize_t>(grid->ghost_layers());
    const exanb::Vec3d sentinel_pos{origin.x + (off.i + gl + 0.5) * cs, origin.y + (off.j + gl + 0.5) * cs,
                                    origin.z + (off.k + gl + 0.5) * cs};
    exanb_naive::insert_sphere(*candidate_grid, sentinel_pos, exanb_naive::sentinel_id, 0.0, 0,
                               std::numeric_limits<uint64_t>::max(), /*confirmed=*/1);

    const size_t n = candidates->size();
    for (size_t i = 0; i < n; i++) {
      const uint64_t id = (uint64_t(rank) << 40) | uint64_t(i);
      exanb_naive::insert_sphere(*candidate_grid, candidates->pos[i], id, candidates->radius[i], candidates->phase[i],
                                 candidates->priority[i]);
    }
    candidate_grid->rebuild_particle_offsets();
  }
};

template <class GridT>
using RSAMPIExanbCandidatesToGridTmpl = RSAMPIExanbCandidatesToGrid<GridT>;

ONIKA_AUTORUN_INIT(exanb_candidates_to_grid) {
  OperatorNodeFactory::instance()->register_factory("exanb_candidates_to_grid",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbCandidatesToGridTmpl>);
}
}  // namespace rsa_mpi
