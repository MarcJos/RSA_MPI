#include <mpi.h>
#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

// RSA MPI
#include <list_of_voxels.hxx>
#include <radius_generator.hxx>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbVoxelAlgorithmContinue : public OperatorNode {
  ADD_SLOT(voxel_list::list_of_voxels<3>, uncovered_voxels, INPUT, REQUIRED,
           DocString{"Voxels not yet covered by any placed sphere"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius generator state"});
  ADD_SLOT(bool, continue_loop, INPUT_OUTPUT, true, DocString{"True as long as another draw should be attempted"});

  inline void execute() override final {
    const int64_t local_nb_vox = int64_t(uncovered_voxels->size());
    int64_t total_nb_vox = 0;
    MPI_Allreduce(&local_nb_vox, &total_nb_vox, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    *continue_loop = (total_nb_vox != 0) && RSARadiusGenerator->is_there_still_radii();
  }
};

template <class GridT>
using RSAMPIExanbVoxelAlgorithmContinueTmpl = RSAMPIExanbVoxelAlgorithmContinue<GridT>;

ONIKA_AUTORUN_INIT(exanb_voxel_algorithm_continue) {
  OperatorNodeFactory::instance()->register_factory(
      "exanb_voxel_algorithm_continue", exanb::make_grid_variant_operator<RSAMPIExanbVoxelAlgorithmContinueTmpl>);
}
}  // namespace rsa_mpi
