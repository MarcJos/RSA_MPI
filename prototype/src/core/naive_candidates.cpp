#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>
#include <rsa_domain.hxx>
#include <vector>

// RSA MPI
#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;
using namespace sac_de_billes;

template <int DIM>
class RSAMPINaiveCandidates : public OperatorNode {
  struct RandomGeneratorState {
    bool initialized = false;
    std::mt19937 rng;
  };

  ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
  ADD_SLOT(int, seed, INPUT, 0);
  ADD_SLOT(int, size, INPUT, REQUIRED, DocString{"Number of shots drawn per MPI process, for a single draw"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});

  ADD_SLOT(rsa_data_storage<DIM>, candidates, OUTPUT, DocString{"Candidate spheres drawn for this iteration"});
  ADD_SLOT(int, shots, OUTPUT, DocString{"Number of shots drawn (locally) that produced these candidates"});

  ADD_SLOT(RandomGeneratorState, random_generator_state, PRIVATE);

  inline std::string documentation() const override final {
    return R"EOF(
Draws candidate sphere centers uniformly over the whole domain, together with
their radius and priority. Meant to be followed by commit_candidates, which
actually adds them to the grid (shared by every drawing strategy).

      YAML example:

        naive_algorithm:
          loop: true
          condition: continue_loop
          body:
            - naive_candidates
            - commit_candidates
            - naive_draw_continue
      )EOF";
  }

  inline void execute() override final {
    rsa_domain<DIM>& domain = *RSADomain;

    if (!random_generator_state->initialized) {
      random_generator_state->rng.seed(*seed);
      random_generator_state->initialized = true;
    }
    std::mt19937& random_generator = random_generator_state->rng;
    auto ulaw = algorithm::auxi::create_random_law(domain, random_generator);
    // get data storage
    auto center_generator = [&ulaw](int size) {
      return algorithm::auxi::generate_sphere_positions<DIM>(ulaw.data(), size);
    };
    //!
    auto priority_generator = [&random_generator](int size) {
      return algorithm::generate_priority<int>(size, random_generator);
    };
    //!
    const auto& radius_gen = *RSARadiusGenerator;
    auto radius_generator = [&radius_gen, &random_generator](int size) { return radius_gen(size, random_generator); };

    *shots = *size;
    *candidates = algorithm::generate_candidates<DIM>(center_generator, radius_generator, priority_generator, *shots);
  }
};

template <int DIM>
using RSAMPINaiveCandidatesTmpl = RSAMPINaiveCandidates<DIM>;

ONIKA_AUTORUN_INIT(naive_candidates) {
  OperatorNodeFactory::instance()->register_factory("naive_candidates",
                                                    make_rsa_mpi_operator<RSAMPINaiveCandidatesTmpl>);
}
}  // namespace rsa_mpi
