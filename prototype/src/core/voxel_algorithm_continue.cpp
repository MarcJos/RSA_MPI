#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <rsa_domain.hxx>
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

template <int DIM>
class RSAMPIVoxelAlgorithmContinue : public OperatorNode {
  ADD_SLOT(voxel_list::list_of_voxels<DIM>, uncovered_voxels, INPUT, REQUIRED,
           DocString{"Voxels not yet covered by any placed sphere"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});
  ADD_SLOT(bool, continue_loop, INPUT_OUTPUT, true, DocString{"True as long as another draw should be attempted"});

  inline std::string documentation() const override final {
    return R"EOF(
Loop condition for the voxel RSA algorithm. Stops when the domain is fully
packed (no uncovered voxel left), or when the radius generator has no more
radii to place.

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
    const int64_t total_nb_vox = rsa_mpi::compute_mpi_sum(int64_t(uncovered_voxels->size()));
    *continue_loop = (total_nb_vox != 0) && RSARadiusGenerator->is_there_still_radii();
  }
};

template <int DIM>
using RSAMPIVoxelAlgorithmContinueTmpl = RSAMPIVoxelAlgorithmContinue<DIM>;

ONIKA_AUTORUN_INIT(voxel_algorithm_continue) {
  OperatorNodeFactory::instance()->register_factory("voxel_algorithm_continue",
                                                    make_rsa_mpi_operator<RSAMPIVoxelAlgorithmContinueTmpl>);
}
}  // namespace rsa_mpi
