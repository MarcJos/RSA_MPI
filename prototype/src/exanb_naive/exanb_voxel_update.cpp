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
#include <RSAMPI/exanb_naive_voxel.hxx>
#include <list_of_voxels.hxx>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbVoxelUpdate : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"exaNBody grid that contains particles"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT, REQUIRED, DocString{"Radius generator state"});
  ADD_SLOT(voxel_list::list_of_voxels<3>, uncovered_voxels, INPUT_OUTPUT,
           DocString{"Voxels not yet covered by any placed sphere"});
  ADD_SLOT(int, shots, INPUT, REQUIRED, DocString{"Number of shots drawn (locally) this iteration"});
  ADD_SLOT(uint64_t, nb_added_spheres_local, INPUT, REQUIRED,
           DocString{"Number of spheres committed (locally) this iteration"});
  ADD_SLOT(bool, verbose, INPUT, false, DocString{"Log the resulting (local max / global total) voxel count"});

  inline void execute() override final {
    const int64_t local_nb_miss = int64_t(*shots) - int64_t(*nb_added_spheres_local);
    int64_t total_nb_miss = 0, total_nb_shots = 0;
    MPI_Allreduce(&local_nb_miss, &total_nb_miss, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    const int64_t local_nb_shots = *shots;
    MPI_Allreduce(&local_nb_shots, &total_nb_shots, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    const double miss_rate = double(total_nb_miss) / (1e-6 + double(total_nb_shots));

    const double desired_miss_rate = algorithm::auxi::magical_default_miss_rate<3>();
    exanb_naive::update_covered_voxels(*uncovered_voxels, *grid, RSARadiusGenerator->get_min_radius(), miss_rate,
                                       desired_miss_rate);

    if (*verbose) {
      int64_t nb_voxels = int64_t(uncovered_voxels->size());
      int64_t max_nb_voxels = 0, total_nb_voxels = 0;
      MPI_Allreduce(&nb_voxels, &max_nb_voxels, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&nb_voxels, &total_nb_voxels, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
      onika::lout << "Total nb voxels: " << total_nb_voxels << " ; Max nb voxels: " << max_nb_voxels << std::endl;
    }
  }
};

template <class GridT>
using RSAMPIExanbVoxelUpdateTmpl = RSAMPIExanbVoxelUpdate<GridT>;

ONIKA_AUTORUN_INIT(exanb_voxel_update) {
  OperatorNodeFactory::instance()->register_factory("exanb_voxel_update",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbVoxelUpdateTmpl>);
}
}  // namespace rsa_mpi
