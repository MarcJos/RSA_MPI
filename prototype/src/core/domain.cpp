#include <vector>
#include <array>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <onika/parallel/parallel_execution_context.h>

#include <rsa_data_storage.hxx>
#include <rsa_random.hxx>
#include <rsa_decoration.hxx>
#include <rsa_domain.hxx>
#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi
{
  using namespace onika;
  using namespace onika::scg;

	template <int DIM> class 
		RSAMPIDomain : public OperatorNode
	{
    using VecND = std::array<double, DIM>;
		ADD_SLOT(RSADim<DIM>, rsa_mpi_dim, INPUT, REQUIRED, DocString{"Used to define the spatial dimension."});
		ADD_SLOT(rsa_domain<DIM>, RSADomain, OUTPUT, DocString{""});
		ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{""}); 
		ADD_SLOT(double, radius, INPUT, REQUIRED, DocString{"Radius of particles"});
		ADD_SLOT(VecND, inf, INPUT, REQUIRED, DocString{"Minimum coordinates of the system. Example inf: [0, ..., 0]"});
		ADD_SLOT(VecND, sup, INPUT, REQUIRED, DocString{"Maximum coordinates of the system, example sup: [0, ..., 0]"});
//		ADD_SLOT(std::array<double, DIM>, inf, INPUT, REQUIRED, DocString{"Minimum coordinates of the system. Example inf: [0, ..., 0]"});
//		ADD_SLOT(std::array<double, DIM>, sup, INPUT, REQUIRED, DocString{"Maximum coordinates of the system, example sup: [0, ..., 0]"});

		inline std::string documentation() const override final
		{
			return R"EOF()EOF";
		}

		inline void execute() override final
		{
			rsa_domain<DIM>& domain = *RSADomain;
			domain = rsa_domain<DIM>(*inf, *sup, *ghost_layer, *radius);
		}
	};

	template<int DIM> using RSAMPIDomainTmpl = RSAMPIDomain<DIM>;

	ONIKA_AUTORUN_INIT(parameters) 
	{ 
		OperatorNodeFactory::instance()->register_factory("rsa_mpi_domain", make_rsa_mpi_operator<RSAMPIDomainTmpl>); 
	} 
}
