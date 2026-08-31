#include <onika/log.h>
#include <onika/parallel/parallel_execution_context.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <rsa_data_storage.hxx>
#include <rsa_decoration.hxx>
#include <rsa_domain.hxx>
#include <rsa_grid.hxx>
#include <rsa_random.hxx>
#include <vector>

// RSA MPI
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <int DIM>
class RSAMPIDisplayGrid : public OperatorNode {
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, REQUIRED, DocString{"Grid that contains particles"});

  inline std::string documentation() const override final { return R"EOF()EOF"; }

  inline void execute() override final {
    const rsa_grid<DIM>& grid = *RSAGrid;
    lout << " ================================== " << std::endl;
    lout << " Number of Spheres  = " << grid.get_number_of_spheres() << std::endl;
    lout << " Spheres Volume     = " << grid.local_volume_of_spheres() << std::endl;
    lout << " ================================== " << std::endl;
  }
};

template <int DIM>
using RSAMPIDisplayGridTmpl = RSAMPIDisplayGrid<DIM>;

ONIKA_AUTORUN_INIT(rsa_mpi_grid) {
  OperatorNodeFactory::instance()->register_factory("rsa_mpi_display_grid",
                                                    make_rsa_mpi_operator<RSAMPIDisplayGridTmpl>);
}
}  // namespace rsa_mpi
