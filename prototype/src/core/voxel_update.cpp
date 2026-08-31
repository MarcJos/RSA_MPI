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
#include <rsa_grid.hxx>
#include <vector>

// RSA MPI
#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;
using namespace sac_de_billes;

template <int DIM>
class RSAMPIVoxelUpdate : public OperatorNode {
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, DocString{"Grid that contains particles"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});
  ADD_SLOT(voxel_list::list_of_voxels<DIM>, uncovered_voxels, INPUT_OUTPUT,
           DocString{"Voxels not yet covered by any placed sphere"});
  ADD_SLOT(int, shots, INPUT, REQUIRED, DocString{"Number of shots drawn (locally) this iteration"});
  ADD_SLOT(uint64_t, nb_added_spheres, INPUT, REQUIRED,
           DocString{"Number of spheres committed (locally) this iteration"});
  ADD_SLOT(bool, verbose, INPUT, false, DocString{"Log the resulting (local max / global total) voxel count"});

  inline std::string documentation() const override final {
    return R"EOF(
Updates the uncovered_voxels list after a draw: removes voxels fully covered by
placed spheres and subdivides the remaining ones, but only once the observed
miss rate (shots that did not turn into a validated sphere) exceeds a default,
dimension-dependent threshold - so that voxels aren't subdivided needlessly.

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
    const int64_t local_nb_miss = int64_t(*shots) - int64_t(*nb_added_spheres);
    const int64_t total_nb_miss = rsa_mpi::compute_mpi_sum(local_nb_miss);
    const int64_t total_nb_shots = rsa_mpi::compute_mpi_sum(int64_t(*shots));
    const double miss_rate = double(total_nb_miss) / (1e-6 + double(total_nb_shots));

    const double desired_miss_rate = algorithm::auxi::magical_default_miss_rate<DIM>();
    algorithm::auxi::update_covered_voxels<DIM>(*uncovered_voxels, *RSAGrid, RSARadiusGenerator->get_min_radius(),
                                                miss_rate, desired_miss_rate, *verbose);
  }
};

template <int DIM>
using RSAMPIVoxelUpdateTmpl = RSAMPIVoxelUpdate<DIM>;

ONIKA_AUTORUN_INIT(voxel_update) {
  OperatorNodeFactory::instance()->register_factory("voxel_update", make_rsa_mpi_operator<RSAMPIVoxelUpdateTmpl>);
}
}  // namespace rsa_mpi
