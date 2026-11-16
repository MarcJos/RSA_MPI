#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <operator_algorithm.hxx>
//
#include <list_of_voxels.hxx>
#include <radius_generator.hxx>
#include <rsa_decoration.hxx>
#include <vector>

// RSA MPI
#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;
using namespace sac_de_billes;

template <int DIM>
class RSAMPIVoxelCandidates : public OperatorNode {
  struct RandomGeneratorState {
    bool initialized = false;
    std::mt19937 rng;
  };

  ADD_SLOT(voxel_list::list_of_voxels<DIM>, uncovered_voxels, INPUT, REQUIRED,
           DocString{"Voxels not yet covered by any placed sphere"});
  ADD_SLOT(int, seed, INPUT, 0);
  ADD_SLOT(int, size, INPUT, REQUIRED,
           DocString{"Ceiling on the (expected) number of shots per draw, across all MPI processes"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});

  ADD_SLOT(rsa_data_storage<DIM>, candidates, OUTPUT, DocString{"Candidate spheres drawn for this iteration"});
  ADD_SLOT(int, shots, OUTPUT, DocString{"Number of shots drawn (locally, Poisson-distributed) this iteration"});

  ADD_SLOT(RandomGeneratorState, random_generator_state, PRIVATE);

  inline std::string documentation() const override final {
    return R"EOF(
Draws candidate sphere centers from the still-uncovered voxels, together with
their radius and priority. Meant to be followed by commit_candidates (adds
them to the grid) and voxel_update (updates the voxel list accordingly).

      YAML example:

        voxel_algorithm:
          loop: true
          condition: continue_loop
          body:
            - voxel_candidates
            - commit_candidates
            - voxel_update
            - voxel_algorithm_continue
      )EOF";
  }

  inline void execute() override final {
    if (!random_generator_state->initialized) {
      random_generator_state->rng.seed(*seed);
      random_generator_state->initialized = true;
    }
    std::mt19937& random_generator = random_generator_state->rng;

    // warning : assumes all the rsa_domains are all equally subdivided with the same voxel_lengths
    const double intensity = algorithm::auxi::compute_intensity_poisson(int64_t(uncovered_voxels->size()));
    const double maximum_intensity = rsa_mpi::compute_mpi_max(intensity);
    *shots = algorithm::auxi::compute_nb_shots_voxel(random_generator, intensity, maximum_intensity, *size);

    const auto& voxels = *uncovered_voxels;
    auto center_generator = [&voxels, &random_generator](int n) { return voxels.pick_points(n, random_generator); };
    //!
    auto priority_generator = [&random_generator](int n) {
      return algorithm::generate_priority<int>(n, random_generator);
    };
    //!
    const auto& radius_gen = *RSARadiusGenerator;
    auto radius_generator = [&radius_gen, &random_generator](int n) { return radius_gen(n, random_generator); };

    *candidates = algorithm::generate_candidates<DIM>(center_generator, radius_generator, priority_generator, *shots);
  }
};

template <int DIM>
using RSAMPIVoxelCandidatesTmpl = RSAMPIVoxelCandidates<DIM>;

ONIKA_AUTORUN_INIT(voxel_candidates) {
  OperatorNodeFactory::instance()->register_factory("voxel_candidates",
                                                    make_rsa_mpi_operator<RSAMPIVoxelCandidatesTmpl>);
}
}  // namespace rsa_mpi
