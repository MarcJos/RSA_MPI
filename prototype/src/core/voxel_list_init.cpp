#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <array>
#include <rsa_domain.hxx>
//
#include <list_of_voxels.hxx>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>
#include <vector>

// RSA MPI
#include <RSAMPI/RSADim.hpp>
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;
using namespace sac_de_billes;

template <int DIM>
class RSAMPIVoxelListInit : public OperatorNode {
  ADD_SLOT(rsa_domain<DIM>, RSADomain, INPUT, REQUIRED, DocString{"RSAMPI Domain"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<DIM>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius of particles"});

  ADD_SLOT(voxel_list::list_of_voxels<DIM>, uncovered_voxels, OUTPUT,
           DocString{"Voxels not yet covered by any placed sphere"});

  inline std::string documentation() const override final {
    return R"EOF(
Builds the initial list of voxels covering the whole domain, sized so that no
voxel's diagonal exceeds max_radius + min_radius. Meant to be run once, before
the voxel_algorithm loop.

      YAML example:

        initialize_algo:
          - rsa_mpi_grid
          - voxel_list_init
      )EOF";
  }

  inline void execute() override final {
    rsa_domain<DIM>& domain = *RSADomain;
    const double max_diagonal = RSARadiusGenerator->get_max_radius() + RSARadiusGenerator->get_min_radius();
    *uncovered_voxels =
        voxel_list::list_of_voxels<DIM>(domain.get_inf(), domain.get_sup() - domain.get_inf(), max_diagonal);
  }
};

template <int DIM>
using RSAMPIVoxelListInitTmpl = RSAMPIVoxelListInit<DIM>;

ONIKA_AUTORUN_INIT(voxel_list_init) {
  OperatorNodeFactory::instance()->register_factory("voxel_list_init", make_rsa_mpi_operator<RSAMPIVoxelListInitTmpl>);
}
}  // namespace rsa_mpi
