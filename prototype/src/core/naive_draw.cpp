#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>
#include <rsa_decoration.hxx>
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
class RSAMPINaiveDraw : public OperatorNode {
  // keeps the random engine alive across successive draws (naive_draw is meant to be
  // called repeatedly from a loop), seeding it only once from *seed
  struct RandomGeneratorState {
    bool initialized = false;
    std::mt19937 rng;
  };

  ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, DocString{"Grid that contains particles"});
  ADD_SLOT(int, seed, INPUT, 0);
  ADD_SLOT(int, size, INPUT, REQUIRED, DocString{"Number of shots drawn per MPI process, for a single draw"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT_OUTPUT, DocString{"Radius of particles"});

  ADD_SLOT(rsa_data_storage<DIM>, ghost_data, INPUT_OUTPUT, DocString{"Store ghost spheres."});
  ADD_SLOT(uint64_t, nb_added_spheres, OUTPUT, DocString{"Number of generated spheres"});

  ADD_SLOT(RandomGeneratorState, random_generator_state, PRIVATE);

  inline std::string documentation() const override final {
    return R"EOF(
 
      YAML example:

        - naive_draw

      )EOF";
  }

  inline void execute() override final {
    rsa_domain<DIM>& domain = *RSADomain;
    rsa_grid<DIM>& grid = *RSAGrid;

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

    // recomputed at each draw, since the radius generator state evolves as spheres get placed
    const uint64_t nb_spheres_total_max = RSARadiusGenerator->get_current_number();
    const int64_t total_nb_shots = int64_t(*size) * rsa_mpi::get_number_of_mpi_processes();
    const bool may_outreach_nb_spheres = (total_nb_shots >= int64_t(nb_spheres_total_max));

    *nb_added_spheres = algorithm::generate_spheres<DIM>(
        grid, domain.get_recv_buffers(), domain.get_send_buffers(), domain.get_ghost_areas(), center_generator,
        radius_generator, priority_generator, *ghost_data, *size, nb_spheres_total_max, may_outreach_nb_spheres);

    // update the radius generator with the (globally) placed spheres, so that
    // the next draw (or the loop continuation check) sees an up to date state
    const uint64_t nb_added_spheres_glob = rsa_mpi::compute_mpi_sum(*nb_added_spheres);
    RSARadiusGenerator->update_placed(nb_added_spheres_glob);
  }
};

template <int DIM>
using RSAMPINaiveDrawTmpl = RSAMPINaiveDraw<DIM>;

ONIKA_AUTORUN_INIT(draw) {
  OperatorNodeFactory::instance()->register_factory("naive_draw", make_rsa_mpi_operator<RSAMPINaiveDrawTmpl>);
}
}  // namespace rsa_mpi
