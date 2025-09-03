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
#include <rsa_grid.hxx>
#include <RSAMPI/make_variant_operator.hpp>

#include <array>

namespace rsa_mpi
{
  using namespace onika;
  using namespace onika::scg;

	template <int DIM> class 
		RSAMPIGrid : public OperatorNode
	{
    using VecND = std::array<double, DIM>;
		ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{""});
		ADD_SLOT(rsa_grid<DIM>, RSAGrid, OUTPUT, DocString{"Grid that contains particles"});

		inline std::string documentation() const override final
		{
			return R"EOF()EOF";
		}

		inline void execute() override final
		{
			rsa_domain<DIM>& domain = *RSADomain;
			rsa_grid<DIM>& grid = *RSAGrid;
			grid = rsa_grid<DIM>(domain.get_m_rad(), domain.get_ghost_layer(), domain.get_inf(), domain.get_sup());
		}
	};


	template<int DIM> using RSAMPIGridTmpl = RSAMPIGrid<DIM>;

	ONIKA_AUTORUN_INIT(rsa_mpi_grid) 
	{ 
		OperatorNodeFactory::instance()->register_factory("rsa_mpi_grid", make_rsa_mpi_operator<RSAMPIGridTmpl>); 
	} 
}
