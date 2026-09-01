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
#include <RSAMPI/fields.h>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbIterationLog : public OperatorNode {
  ADD_SLOT(exanb::Domain, domain, INPUT, REQUIRED, DocString{"exaNBody domain (global bounds)"});
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"exaNBody grid"});
  ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{"Thickness (in cells) of the grid's ghost shell"});
  ADD_SLOT(uint64_t, nb_added_spheres, INPUT, REQUIRED, DocString{"Global number of spheres committed this round"});

  ADD_SLOT(int, iteration, PRIVATE, 0);

  inline void execute() override final {
    ++(*iteration);

    const ssize_t gl = *ghost_layer;
    const exanb::IJK dim = grid->dimension();
    uint64_t local_particles = 0;
    double local_volume = 0.0;
    for (ssize_t i = gl; i < dim.i - gl; i++) {
      for (ssize_t j = gl; j < dim.j - gl; j++) {
        for (ssize_t k = gl; k < dim.k - gl; k++) {
          auto& cell = grid->cell(exanb::IJK{i, j, k});
          const size_t n = cell.size();
          local_particles += n;
          for (size_t s = 0; s < n; s++) {
            const double r = cell[exanb::field::radius][s];
            local_volume += (4.0 / 3.0) * M_PI * r * r * r;
          }
        }
      }
    }

    uint64_t nb_particles = 0;
    double total_spheres_volume = 0.0;
    MPI_Allreduce(&local_particles, &nb_particles, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_volume, &total_spheres_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const exanb::Vec3d extent = domain->bounds_size();
    const double total_volume = extent.x * extent.y * extent.z;
    const double compacity = total_spheres_volume / total_volume;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      onika::lout << " [iteration " << *iteration << "] particles = " << nb_particles
                  << " ; added = " << *nb_added_spheres << " ; compacity = " << compacity << std::endl;
    }
  }
};

template <class GridT>
using RSAMPIExanbIterationLogTmpl = RSAMPIExanbIterationLog<GridT>;

ONIKA_AUTORUN_INIT(exanb_iteration_log) {
  OperatorNodeFactory::instance()->register_factory("exanb_iteration_log",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbIterationLogTmpl>);
}
}  // namespace rsa_mpi
