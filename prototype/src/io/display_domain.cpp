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
#include <rsa_parameters.hxx>
#include <rsa_domain.hxx>
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi
{
  using namespace onika;
  using namespace onika::scg;

	template <int DIM> class RSAMPIDisplayDomain : public OperatorNode
	{
    ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{""});

		inline std::string documentation() const override final
		{
			return R"EOF()EOF";
		}

		inline void execute() override final
		{
      rsa_domain<DIM>& domain = *RSADomain;
      domain.domain_log();
		}
	};


  template<int DIM> using RSAMPIDisplayDomainTmpl = RSAMPIDisplayDomain<DIM>;

  ONIKA_AUTORUN_INIT(parameters) 
  { 
    OperatorNodeFactory::instance()->register_factory("rsa_mpi_display_domain", make_rsa_mpi_operator<RSAMPIDisplayDomainTmpl>); 
  } 
}
