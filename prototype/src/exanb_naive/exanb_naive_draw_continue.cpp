#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

// RSA MPI
#include <radius_generator.hxx>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbNaiveDrawContinue : public OperatorNode {
  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius generator state"});
  ADD_SLOT(int, n_draw, INPUT, REQUIRED, DocString{"Maximum number of draws for the naive method"});
  ADD_SLOT(int, draw_count, PRIVATE, 0);
  ADD_SLOT(bool, continue_loop, INPUT_OUTPUT, true, DocString{"True as long as another draw should be attempted"});

  inline void execute() override final {
    ++(*draw_count);
    *continue_loop = (*draw_count < *n_draw) && RSARadiusGenerator->is_there_still_radii();
  }
};

template <class GridT>
using RSAMPIExanbNaiveDrawContinueTmpl = RSAMPIExanbNaiveDrawContinue<GridT>;

ONIKA_AUTORUN_INIT(exanb_naive_draw_continue) {
  OperatorNodeFactory::instance()->register_factory(
      "exanb_naive_draw_continue", exanb::make_grid_variant_operator<RSAMPIExanbNaiveDrawContinueTmpl>);
}
}  // namespace rsa_mpi
