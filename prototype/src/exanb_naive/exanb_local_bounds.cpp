#include <onika/log.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <cmath>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbLocalBounds : public OperatorNode {
  ADD_SLOT(exanb::Domain, domain, INPUT, REQUIRED, DocString{"exaNBody domain (global, for total volume)"});
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"exaNBody grid, already partitioned and padded"});
  ADD_SLOT(int, global_size, INPUT, REQUIRED,
           DocString{"Global number of candidates drawn per round, across all ranks"});

  ADD_SLOT(exanb::Vec3d, inf, OUTPUT, DocString{"Minimum coordinates of this rank's real block"});
  ADD_SLOT(exanb::Vec3d, sup, OUTPUT, DocString{"Maximum coordinates of this rank's real block"});
  ADD_SLOT(int, size, OUTPUT, DocString{"Number of candidates this rank should draw this round"});

  inline void execute() override final {
    const ssize_t gl = static_cast<ssize_t>(grid->ghost_layers());
    const exanb::IJK off = grid->offset();
    const exanb::IJK dim = grid->dimension();
    const double cs = grid->cell_size();
    const exanb::Vec3d origin = grid->origin();

    *inf = exanb::Vec3d{origin.x + (off.i + gl) * cs, origin.y + (off.j + gl) * cs, origin.z + (off.k + gl) * cs};
    *sup = exanb::Vec3d{origin.x + (off.i + dim.i - gl) * cs, origin.y + (off.j + dim.j - gl) * cs,
                        origin.z + (off.k + dim.k - gl) * cs};

    const exanb::Vec3d local_extent = *sup - *inf;
    const double local_volume = local_extent.x * local_extent.y * local_extent.z;

    const exanb::Vec3d global_extent = domain->bounds_size();
    const double global_volume = global_extent.x * global_extent.y * global_extent.z;

    *size = static_cast<int>(std::lround((*global_size) * (local_volume / global_volume)));
    if (*size < 0) *size = 0;
  }
};

template <class GridT>
using RSAMPIExanbLocalBoundsTmpl = RSAMPIExanbLocalBounds<GridT>;

ONIKA_AUTORUN_INIT(exanb_local_bounds) {
  OperatorNodeFactory::instance()->register_factory("exanb_local_bounds",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbLocalBoundsTmpl>);
}
}  // namespace rsa_mpi
