#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>
#include <array>
#include <radius_generator.hxx>
#include <rsa_domain.hxx>
#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

// Generalizes uniform_radius_generator to several sphere sizes drawn in the same
// run (e.g. a bimodal packing): each radius targets its own volume fraction of
// the total domain volume, largest radius first.
//
// Does not take the domain itself (only its inf/sup) so that it can run before
// rsa_mpi_domain: rsa_mpi_domain's own cell_size can then be wired directly to
// this operator's cell_size output, instead of duplicating the largest radius.
template <int DIM>
class RSAMPIMultiRadiusGenerator : public OperatorNode {
  using VecND = std::array<double, DIM>;
  // anchors this operator's DIM (no RSADomain/RSADim-typed input would otherwise
  // pin it down, since this operator is meant to run before rsa_mpi_domain)
  ADD_SLOT(RSADim<DIM>, rsa_mpi_dim, INPUT, REQUIRED, DocString{"Used to define the spatial dimension."});
  ADD_SLOT(VecND, inf, INPUT, REQUIRED, DocString{"Minimum coordinates of the system. Example inf: [0, ..., 0]"});
  ADD_SLOT(VecND, sup, INPUT, REQUIRED, DocString{"Maximum coordinates of the system, example sup: [0, ..., 0]"});
  ADD_SLOT(std::vector<double>, radii, INPUT, REQUIRED,
           DocString{"Sphere radii, strictly decreasing (largest first)"});
  ADD_SLOT(std::vector<double>, volume_fractions, INPUT, REQUIRED,
           DocString{"Desired volume fraction of the total domain volume for each radius, in ]0,1]"});
  ADD_SLOT(std::vector<int>, phases, INPUT, std::vector<int>{},
           DocString{"Phase id for each radius (defaults to 0 for all)"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, OUTPUT, DocString{"Radius of particles"});
  ADD_SLOT(double, cell_size, OUTPUT,
           DocString{"Largest configured radius; meant to feed rsa_mpi_domain's cell_size, so the domain's "
                     "implicit cell size always covers the largest sphere without duplicating the value"});

  inline std::string documentation() const override final {
    return R"EOF(
Builds a radius generator drawing several sphere sizes in the same run, each
targeting its own volume fraction of the total domain volume. Also exposes
the largest configured radius as the "cell_size" output slot, meant to feed
rsa_mpi_domain (positioned after this operator).

      YAML example:

        input_data:
          - multi_radius_generator:
             inf: [0,0,0]
             sup: [1,1,1]
             radii: [0.1, 0.05]
             volume_fractions: [0.25, 1.0]
          - rsa_mpi_domain:
             inf: [0,0,0]
             sup: [1,1,1]
      )EOF";
  }

  inline void execute() override final {
    if (radii->empty()) {
      fatal_error() << "multi_radius_generator: radii must not be empty" << std::endl;
    }
    if (radii->size() != volume_fractions->size()) {
      fatal_error() << "multi_radius_generator: radii (" << radii->size() << ") and volume_fractions ("
                    << volume_fractions->size() << ") must have the same size" << std::endl;
    }
    if (!phases->empty() && phases->size() != radii->size()) {
      fatal_error() << "multi_radius_generator: phases (" << phases->size() << ") must be empty or have the same size as radii ("
                    << radii->size() << ")" << std::endl;
    }
    // the underlying RadiusGenerator takes the *first* entry as its max_radius and the
    // *last* one as its min_radius (see radius_generator.ixx); if radii isn't sorted in
    // strictly decreasing order, that assumption silently breaks and the cell_size
    // deduced below would not actually cover the largest configured sphere.
    for (size_t i = 1; i < radii->size(); i++) {
      if ((*radii)[i] >= (*radii)[i - 1]) {
        fatal_error() << "multi_radius_generator: radii must be strictly decreasing (largest first); got "
                      << (*radii)[i - 1] << " followed by " << (*radii)[i] << std::endl;
      }
    }

    const double volume = sac_de_billes::auxi_function::productOf<double>(*sup - *inf);

    std::vector<std::tuple<double, double, int>> desired_radius_volumeFraction_phase(radii->size());
    for (size_t i = 0; i < radii->size(); i++) {
      const int phase = phases->empty() ? 0 : (*phases)[i];
      desired_radius_volumeFraction_phase[i] = {(*radii)[i], (*volume_fractions)[i], phase};
    }
    *RSARadiusGenerator = sac_de_billes::RadiusGenerator<DIM>(desired_radius_volumeFraction_phase, volume);
    *cell_size = radii->front();

    onika::lout << " ================================== " << std::endl;
    onika::lout << " Radius Generator Mode = Multi " << std::endl;
    for (size_t i = 0; i < radii->size(); i++) {
      onika::lout << " Radius[" << i << "] = " << (*radii)[i] << " ; Volume fraction = " << (*volume_fractions)[i]
                  << std::endl;
    }
    onika::lout << " Volume = " << volume << std::endl;
    onika::lout << " ================================== " << std::endl;
  }
};

template <int DIM>
using RSAMPIMultiRadiusGeneratorTmpl = RSAMPIMultiRadiusGenerator<DIM>;

ONIKA_AUTORUN_INIT(multi_radius_generator) {
  OperatorNodeFactory::instance()->register_factory("multi_radius_generator",
                                                     make_rsa_mpi_operator<RSAMPIMultiRadiusGeneratorTmpl>);
}
}  // namespace rsa_mpi
