#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

// RSA MPI
#include <list_of_voxels.hxx>
#include <radius_generator.hxx>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbVoxelListInit : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"exaNBody grid, already partitioned and padded"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius generator state"});

  ADD_SLOT(voxel_list::list_of_voxels<3>, uncovered_voxels, OUTPUT,
           DocString{"Voxels not yet covered by any placed sphere, over this rank's local block"});

  inline void execute() override final {
    const ssize_t gl = static_cast<ssize_t>(grid->ghost_layers());
    const exanb::IJK off = grid->offset();
    const exanb::IJK dim = grid->dimension();
    const double cs = grid->cell_size();
    const exanb::Vec3d origin = grid->origin();

    const exanb::Vec3d bmin{origin.x + (off.i + gl) * cs, origin.y + (off.j + gl) * cs, origin.z + (off.k + gl) * cs};
    const exanb::Vec3d bmax{origin.x + (off.i + dim.i - gl) * cs, origin.y + (off.j + dim.j - gl) * cs,
                            origin.z + (off.k + dim.k - gl) * cs};

    const sac_de_billes::Point<3> inf{bmin.x, bmin.y, bmin.z};
    const sac_de_billes::Point<3> extent{bmax.x - bmin.x, bmax.y - bmin.y, bmax.z - bmin.z};
    const double max_diagonal = RSARadiusGenerator->get_max_radius() + RSARadiusGenerator->get_min_radius();
    *uncovered_voxels = voxel_list::list_of_voxels<3>(inf, extent, max_diagonal);
  }
};

template <class GridT>
using RSAMPIExanbVoxelListInitTmpl = RSAMPIExanbVoxelListInit<GridT>;

ONIKA_AUTORUN_INIT(exanb_voxel_list_init) {
  OperatorNodeFactory::instance()->register_factory("exanb_voxel_list_init",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbVoxelListInitTmpl>);
}
}  // namespace rsa_mpi
