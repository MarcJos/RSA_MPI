#include <onika/log.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <array>
#include <cmath>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbDomainInit : public OperatorNode {
  using BoolND = std::array<bool, 3>;

  ADD_SLOT(exanb::Vec3d, domain_inf, INPUT, REQUIRED,
           DocString{"Minimum coordinates of the system. Example domain_inf: [0, 0, 0]"});
  ADD_SLOT(exanb::Vec3d, domain_sup, INPUT, REQUIRED,
           DocString{"Maximum coordinates of the system, example domain_sup: [1, 1, 1]"});
  ADD_SLOT(double, cell_size, INPUT, REQUIRED, DocString{"Implicit cell size (must be >= the largest sphere radius)"});
  ADD_SLOT(BoolND, periodic, INPUT, BoolND{true, true, true},
           DocString{"Per-axis periodicity. Example periodic: [true, true, false]"});

  ADD_SLOT(exanb::Domain, domain, INPUT_OUTPUT,
           DocString{"exaNBody domain (the global, unpartitioned [inf,sup] region)"});
  ADD_SLOT(GridT, grid, OUTPUT, DocString{"empty exaNBody grid - init_rcb_grid requires an empty grid to partition"});

  inline void execute() override final {
    const exanb::Vec3d& bmin = *domain_inf;
    const exanb::Vec3d& bmax = *domain_sup;

    const double actual_cell_size = 2.0 * (*cell_size);

    domain->set_bounds(exanb::AABB{bmin, bmax});
    domain->set_cell_size(actual_cell_size);
    domain->set_periodic_boundary((*periodic)[0], (*periodic)[1], (*periodic)[2]);

    const exanb::Vec3d extent = bmax - bmin;
    const exanb::IJK real_dims{static_cast<ssize_t>(std::ceil(extent.x / actual_cell_size)),
                               static_cast<ssize_t>(std::ceil(extent.y / actual_cell_size)),
                               static_cast<ssize_t>(std::ceil(extent.z / actual_cell_size))};
    domain->set_grid_dimension(real_dims);
    grid->reset();
  }
};

template <class GridT>
using RSAMPIExanbDomainInitTmpl = RSAMPIExanbDomainInit<GridT>;

ONIKA_AUTORUN_INIT(exanb_domain_init) {
  OperatorNodeFactory::instance()->register_factory("exanb_domain_init",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbDomainInitTmpl>);
}
}  // namespace rsa_mpi
