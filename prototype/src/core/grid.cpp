#include <onika/log.h>
#include <onika/parallel/parallel_execution_context.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <RSAMPI/make_variant_operator.hpp>
#include <array>
#include <radius_generator.hxx>
#include <rsa_data_storage.hxx>
#include <rsa_decoration.hxx>
#include <rsa_domain.hxx>
#include <rsa_grid.hxx>
#include <rsa_random.hxx>
#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <int DIM>
class RSAMPIGrid : public OperatorNode {
  using VecND = std::array<double, DIM>;
  ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{""});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, OUTPUT, DocString{"Grid that contains particles"});

  inline std::string documentation() const override final { return R"EOF()EOF"; }

  inline void execute() override final {
    rsa_domain<DIM>& domain = *RSADomain;
    rsa_grid<DIM>& grid = *RSAGrid;
    const double max_radius = RSARadiusGenerator->get_max_radius();
    if (max_radius > domain.get_m_rad()) {
      fatal_error() << "Impossible to have the maximal radius of the radius generator (" << max_radius
                     << ") larger than the domain implicit radius (" << domain.get_m_rad() << ")" << std::endl;
    }
    grid = rsa_grid<DIM>(max_radius, domain.get_ghost_layer(), domain.get_inf(), domain.get_sup());
  }
};

template <int DIM>
using RSAMPIGridTmpl = RSAMPIGrid<DIM>;

ONIKA_AUTORUN_INIT(rsa_mpi_grid) {
  OperatorNodeFactory::instance()->register_factory("rsa_mpi_grid", make_rsa_mpi_operator<RSAMPIGridTmpl>);
}
}  // namespace rsa_mpi
