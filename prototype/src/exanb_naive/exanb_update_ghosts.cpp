#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/mpi/update_ghosts.h>
#include <onika/scg/operator_factory.h>
// RSA MPI
#include <RSAMPI/fields.h>

namespace rsa_mpi {
using RSAGhostFieldSet = ::exanb::FieldSet<::exanb::field::_rx, ::exanb::field::_ry, ::exanb::field::_rz,
                                           ::exanb::field::_id, ::exanb::field::_radius, ::exanb::field::_phase,
                                           ::exanb::field::_priority, ::exanb::field::_confirmed>;

template <class GridT>
using RSAMPIExanbUpdateGhosts = ::exanb::UpdateGhostsNode<GridT, RSAGhostFieldSet, true, false>;

ONIKA_AUTORUN_INIT(exanb_update_ghosts) {
  ::onika::scg::OperatorNodeFactory::instance()->register_factory(
      "exanb_update_ghosts", ::exanb::make_grid_variant_operator<RSAMPIExanbUpdateGhosts>);
}
}  // namespace rsa_mpi
