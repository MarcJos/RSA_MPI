#pragma once

#include<rsa_parameters.hxx>

namespace rsa_mpi
{
	template <int DIM> class RSAMPIParameters : public OperatorNode
	{
    ADD_SLOT(rsa_parameters, RSAParameters, OUTPUT, DocString{""});
    ADD_SLOT(std::array<double, DIM>, inf, INPUT, REQUIRED, DocString{"Minimum coordinates of the system. Example inf: [0, ..., 0]"});
    ADD_SLOT(std::array<double, DIM>, sup, INPUT, REQUIRED, DocString{"aximum coordinates of the system, example sup: [0, ..., 0]"});
    ADD_SLOT(double, radius, INPUT, REQUIRED, DocString{"Radius of particles"});
    ADD_SLOT(bool, paraview, INTPUT, false, DocString{" Flag for Paraview visualization"});
    ADD_SLOT(int, size, INPUT, 6000, DocString{"TO DO RENAME"});
    ADD_SLOT(int, n_draw, INPUT, 10, DocString{"Number of draws"});
    ADD_SLOT(int, seed, INTPUT, 0, DocString{"Seed for random number generation"}); 

		inline std::string documentation() const override final
		{
			return R"EOF()EOF";
		}

		inline void execute() override final
		{
      RSAParameters& params = *RSAParameters;
      params.DIM = DIM;
      params.radius = *radius;
      params.l_min = *inf;
      params.l_max = *sup;
      params.seed = *seed;
      params.n_draw = *n_draw;
      params.size = *size;
      params.paraview = *paraview;      
      params.minimal_requirement();
      params.Display();
		}
	}

  ONIKA_AUTORUN_INIT(parameters) 
  { 
    OperatorNodeFactory::instance()->register_factory("paramters_dim_2", make_simple_operator<RSAMPIParameters<2>>); 
    OperatorNodeFactory::instance()->register_factory("paramters_dim_3", make_simple_operator<RSAMPIParameters<3>>); 
    OperatorNodeFactory::instance()->register_factory("paramters_dim_4", make_simple_operator<RSAMPIParameters<4>>); 
    OperatorNodeFactory::instance()->register_factory("paramters_dim_5", make_simple_operator<RSAMPIParameters<5>>); 
    OperatorNodeFactory::instance()->register_factory("paramters_dim_6", make_simple_operator<RSAMPIParameters<6>>); 
  } 
}
