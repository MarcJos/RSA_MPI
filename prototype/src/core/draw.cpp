#include <vector>
#include <array>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>

#include <radius_generator.hxx>
#include <rsa_domain.hxx>
#include <operator_algorithm.hxx>

#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

#include <array>

namespace rsa_mpi
{
  using namespace onika;
  using namespace onika::scg;
  using namespace sac_de_billes;

	template <int DIM> class 
		RSAMPINaiveDraw : public OperatorNode
	{
		ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
		ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, DocString{"Grid that contains particles"});
    ADD_SLOT(int, seed, INPUT, 0);
    ADD_SLOT(int, nb_shots, INPUT, REQUIRED, DocString{"Number of generated spheres"});
    ADD_SLOT(uint64_t, nb_spheres_total_max, INPUT, REQUIRED, DocString{"Maximum number of generated spheres"});
    ADD_SLOT(uint64_t, may_outreach_nb_spheres, INPUT, REQUIRED, DocString{""});
		ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, DocString{"Radius of particles"});

    ADD_SLOT(rsa_data_storage<DIM>, ghost_data, INPUT_OUTPUT, DocString{"Store ghost spheres."});
    ADD_SLOT(uint64_t, nb_added_spheres, OUTPUT, DocString{"Number of generated spheres"});

		inline std::string documentation() const override final
		{
			return R"EOF(
 
      YAML example:

        - draw

      )EOF";
		}

		inline void execute() override final
		{
			rsa_domain<DIM>& domain = *RSADomain;
  		rsa_grid<DIM>& grid = *RSAGrid;

      // first version
      std::mt19937 random_generator(*seed);
      auto ulaw = algorithm::auxi::create_random_law(domain, random_generator);
      // get data storage
      auto center_generator = [&ulaw](int size) {
        return algorithm::auxi::generate_sphere_positions<DIM>(ulaw.data(), size);
      };
      //!
      auto priority_generator = [&random_generator](int size) {
        return algorithm::generate_priority<int>(size, random_generator);
      };
      //!
      const auto& radius_gen = *RSARadiusGenerator;
      auto radius_generator = [&radius_gen, &random_generator](int size) {
        return radius_gen(size, random_generator);
      };

      *nb_added_spheres = algorithm::generate_spheres<DIM>(
          grid,
          domain.get_recv_buffers(),
          domain.get_send_buffers(),
          domain.get_ghost_areas(),
          center_generator,
          radius_generator,
          priority_generator,          
          *ghost_data,
          *nb_shots,
          *nb_spheres_total_max,
          *may_outreach_nb_spheres);
      // some checks used for debugging
      assert(this->get_grid().check_particles() && " collision detected");
      assert(check_no_doublon(this->get_grid()));
    }
  };


  template<int DIM> using RSAMPINaiveDrawTmpl = RSAMPINaiveDraw<DIM>;

  ONIKA_AUTORUN_INIT(draw) 
  { 
    OperatorNodeFactory::instance()->register_factory("naive_draw", make_rsa_mpi_operator<RSAMPINaiveDrawTmpl>); 
  } 
}
