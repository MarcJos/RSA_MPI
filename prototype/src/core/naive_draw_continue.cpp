#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>
#include <array>
#include <radius_generator.hxx>
#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <int DIM>
class RSAMPINaiveDrawContinue : public OperatorNode {
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});
  ADD_SLOT(int, n_draw, INPUT, REQUIRED, DocString{"Maximum number of draws for the naive method"});
  ADD_SLOT(int, draw_count, PRIVATE, 0);
  // must be INPUT_OUTPUT (not plain OUTPUT): onika's batch "condition:" mechanism only
  // reuses an already-existing same-named slot if one was created during the normal
  // input-wiring pass, which only happens for INPUT/INPUT_OUTPUT slots. A plain OUTPUT
  // here would leave the loop condition permanently disconnected from this value.
  ADD_SLOT(bool, continue_loop, INPUT_OUTPUT, true, DocString{"True as long as another draw should be attempted"});

  inline std::string documentation() const override final {
    return R"EOF(
Loop condition for the naive RSA algorithm. Stops when the maximum number of
draws is reached, or when the radius generator has no more radii to place.

      YAML example:

        naive_algorithm:
          loop: true
          condition: continue_loop
          body:
            - naive_draw
            - naive_draw_continue
      )EOF";
  }

  inline void execute() override final {
    ++(*draw_count);
    *continue_loop = (*draw_count < *n_draw) && RSARadiusGenerator->is_there_still_radii();
  }
};

template <int DIM>
using RSAMPINaiveDrawContinueTmpl = RSAMPINaiveDrawContinue<DIM>;

ONIKA_AUTORUN_INIT(naive_draw_continue) {
  OperatorNodeFactory::instance()->register_factory("naive_draw_continue",
                                                    make_rsa_mpi_operator<RSAMPINaiveDrawContinueTmpl>);
}
}  // namespace rsa_mpi
