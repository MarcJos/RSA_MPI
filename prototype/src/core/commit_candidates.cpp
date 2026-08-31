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
class RSAMPICommitCandidates : public OperatorNode {
  ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, DocString{"Grid that contains particles"});
  ADD_SLOT(rsa_data_storage<DIM>, candidates, INPUT, REQUIRED, DocString{"Candidate spheres to commit"});
  ADD_SLOT(int, size, INPUT, REQUIRED,
           DocString{"Number of shots drawn per MPI process that produced these candidates"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT_OUTPUT, DocString{"Radius of particles"});

  ADD_SLOT(uint64_t, nb_added_spheres, OUTPUT, DocString{"Number of newly committed spheres"});

  inline std::string documentation() const override final {
    return R"EOF(
Commits candidate spheres (produced by a drawing strategy such as
naive_candidates) into the grid, resolving conflicts, then updates the radius
generator with the spheres that got placed.

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
    rsa_grid<DIM>& grid = *RSAGrid;

    // recomputed at each draw, since the radius generator state evolves as spheres get placed
    const uint64_t nb_spheres_total_max = RSARadiusGenerator->get_current_number();
    const int64_t total_nb_shots = int64_t(*size) * rsa_mpi::get_number_of_mpi_processes();
    const bool may_outreach_nb_spheres = (total_nb_shots >= int64_t(nb_spheres_total_max));

    *nb_added_spheres = algorithm::commit_candidates<DIM>(grid, domain.get_recv_buffers(), domain.get_send_buffers(),
                                                          domain.get_ghost_areas(), *candidates, nb_spheres_total_max,
                                                          may_outreach_nb_spheres);

    // update the radius generator with the (globally) placed spheres, so that
    // the next draw (or the loop continuation check) sees an up to date state
    const uint64_t nb_added_spheres_glob = rsa_mpi::compute_mpi_sum(*nb_added_spheres);
    RSARadiusGenerator->update_placed(nb_added_spheres_glob);
  }
};

template <int DIM>
using RSAMPICommitCandidatesTmpl = RSAMPICommitCandidates<DIM>;

ONIKA_AUTORUN_INIT(commit_candidates) {
  OperatorNodeFactory::instance()->register_factory("commit_candidates",
                                                    make_rsa_mpi_operator<RSAMPICommitCandidatesTmpl>);
}
}  // namespace rsa_mpi
