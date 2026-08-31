#include <onika/log.h>
#include <onika/parallel/parallel_execution_context.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <RSAMPI/RSADim.hpp>
#include <array>
#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <int DIM>
class RSAMPISetDim : public OperatorNode {
  ADD_SLOT(RSADim<DIM>, rsa_mpi_dim, OUTPUT, DocString{"Used to define the spatial dimension."});

  inline std::string documentation() const override final { return R"EOF()EOF"; }

  inline void execute() override final {
    lout << "=================================" << std::endl;
    lout << "Dimension: " << DIM << std::endl;
    lout << "=================================" << std::endl;
  }
};

// template<int DIM> using RSAMPISetDimTmpl = RSAMPISetDim<DIM>;

ONIKA_AUTORUN_INIT(set_dim) {
  OperatorNodeFactory::instance()->register_factory("dim_2", make_simple_operator<RSAMPISetDim<2>>);
  OperatorNodeFactory::instance()->register_factory("dim_3", make_simple_operator<RSAMPISetDim<3>>);
  OperatorNodeFactory::instance()->register_factory("dim_4", make_simple_operator<RSAMPISetDim<4>>);
  OperatorNodeFactory::instance()->register_factory("dim_5", make_simple_operator<RSAMPISetDim<5>>);
  OperatorNodeFactory::instance()->register_factory("dim_6", make_simple_operator<RSAMPISetDim<6>>);
  OperatorNodeFactory::instance()->register_factory("dim_7", make_simple_operator<RSAMPISetDim<7>>);
  OperatorNodeFactory::instance()->register_factory("dim_8", make_simple_operator<RSAMPISetDim<8>>);
  OperatorNodeFactory::instance()->register_factory("dim_9", make_simple_operator<RSAMPISetDim<9>>);
  OperatorNodeFactory::instance()->register_factory("dim_10", make_simple_operator<RSAMPISetDim<10>>);
}
}  // namespace rsa_mpi
