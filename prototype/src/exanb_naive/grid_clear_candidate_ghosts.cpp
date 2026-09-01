#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/grid.h>
#include <exanb/core/grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIGridClearCandidateGhosts : public OperatorNode {
  ADD_SLOT(GridT, candidate_grid, INPUT_OUTPUT, DocString{"Grid whose ghost cells get emptied"});

  inline void execute() override final {
    auto& g = *candidate_grid;
    exanb::apply_grid_shell(candidate_grid->dimension(), 0, candidate_grid->ghost_layers(),
                            [&g](ssize_t i, const exanb::IJK&) { g.cell(i).clear(g.cell_allocator()); });
  }
};

template <class GridT>
using RSAMPIGridClearCandidateGhostsTmpl = RSAMPIGridClearCandidateGhosts<GridT>;

ONIKA_AUTORUN_INIT(grid_clear_candidate_ghosts) {
  OperatorNodeFactory::instance()->register_factory(
      "grid_clear_candidate_ghosts", exanb::make_grid_variant_operator<RSAMPIGridClearCandidateGhostsTmpl>);
}
}  // namespace rsa_mpi
