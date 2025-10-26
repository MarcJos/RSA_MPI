#include <vector>
#include <array>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>

#include <radius_generator.hxx>
#include <rsa_domain.hxx>

#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

#include <array>

namespace rsa_mpi
{
  using namespace onika;
  using namespace onika::scg;

	template <int DIM> class 
		RSAMPIUniformRadiusGenerator : public OperatorNode
	{
		ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
		ADD_SLOT(double, radius, INPUT, REQUIRED, DocString{"Radius of particles"});
		ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, OUTPUT, DocString{"Radius of particles"});

		inline std::string documentation() const override final
		{
			return R"EOF()EOF";
		}

		inline void execute() override final
		{
			rsa_domain<DIM>& domain = *RSADomain;
			*RSARadiusGenerator = RadiusGenerator<DIM>(std::vector<std::tuple<double,double, int>>{ {*radius, 1.0, 0} }, domain.get_total_volume());

      onika::lout << " ================================== " << std::endl;
      onika::lout << " Radius Generator Mode = Uniform " << std::endl;
      onika::lout << " Radius = " << *radius << std::endl;
      onika::lout << " Volume = " << domain.get_total_volume() << std::endl;
      onika::lout << " ================================== " << std::endl;
		}
	};


	template<int DIM> using RSAMPIUniformRadiusGeneratorTmpl = RSAMPIUniformRadiusGenerator<DIM>;

	ONIKA_AUTORUN_INIT(uniform_radius_generator) 
	{ 
		OperatorNodeFactory::instance()->register_factory("uniform_radius_generator", make_rsa_mpi_operator<RSAMPIUniformRadiusGeneratorTmpl>); 
	} 
}
