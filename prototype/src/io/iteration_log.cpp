#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <rsa_decoration.hxx>
#include <rsa_domain.hxx>
#include <rsa_grid.hxx>

#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

// Logs, at each algorithm iteration, the running number of placed spheres,
// the number added by this iteration, and the resulting packing fraction
// (compacity = total sphere volume / domain volume).
template <int DIM>
class RSAMPIIterationLog : public OperatorNode {
  ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, DocString{"Grid that contains particles"});
  ADD_SLOT(uint64_t, nb_added_spheres, INPUT, REQUIRED, DocString{"Number of spheres committed (locally) this iteration"});

  ADD_SLOT(int, iteration, PRIVATE, 0);

  inline std::string documentation() const override final {
    return R"EOF(
Logs, once per algorithm iteration, the running number of placed spheres, the
number added this iteration, and the resulting packing fraction (compacity).

      YAML example:

        naive_algorithm:
          loop: true
          condition: continue_loop
          body:
            - naive_candidates
            - commit_candidates
            - rsa_mpi_iteration_log
            - naive_draw_continue
      )EOF";
  }

  inline void execute() override final {
    ++(*iteration);

    const rsa_domain<DIM>& domain = *RSADomain;
    const rsa_grid<DIM>& grid = *RSAGrid;

    const uint64_t nb_added_glob = rsa_mpi::compute_mpi_sum(*nb_added_spheres);
    const uint64_t nb_particles = rsa_mpi::compute_mpi_sum(grid.get_number_of_spheres_fast());
    const double total_spheres_volume = rsa_mpi::compute_mpi_sum(grid.local_volume_of_spheres());
    const double compacity = total_spheres_volume / domain.get_total_volume();

    onika::lout << " [iteration " << *iteration << "] particles = " << nb_particles << " ; added = " << nb_added_glob
                << " ; compacity = " << compacity << std::endl;
  }
};

template <int DIM>
using RSAMPIIterationLogTmpl = RSAMPIIterationLog<DIM>;

ONIKA_AUTORUN_INIT(iteration_log) {
  OperatorNodeFactory::instance()->register_factory("rsa_mpi_iteration_log",
                                                     make_rsa_mpi_operator<RSAMPIIterationLogTmpl>);
}
}  // namespace rsa_mpi
