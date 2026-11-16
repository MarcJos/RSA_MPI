#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbPadGridGhost : public OperatorNode {
  ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{"Thickness (in cells) of the grid's ghost shell"});
  ADD_SLOT(GridT, grid, INPUT_OUTPUT, DocString{"exaNBody grid, freshly partitioned by init_rcb_grid"});

  inline void execute() override final {
    const ssize_t gl = *ghost_layer;
    const exanb::IJK local_offset = grid->offset();
    const exanb::IJK local_dims = grid->dimension();

    grid->set_offset(exanb::IJK{local_offset.i - gl, local_offset.j - gl, local_offset.k - gl});
    grid->set_dimension(exanb::IJK{local_dims.i + 2 * gl, local_dims.j + 2 * gl, local_dims.k + 2 * gl});
    grid->set_max_neighbor_distance(gl * grid->cell_size());
  }
};

template <class GridT>
using RSAMPIExanbPadGridGhostTmpl = RSAMPIExanbPadGridGhost<GridT>;

ONIKA_AUTORUN_INIT(exanb_pad_grid_ghost) {
  OperatorNodeFactory::instance()->register_factory("exanb_pad_grid_ghost",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbPadGridGhostTmpl>);
}
}  // namespace rsa_mpi
