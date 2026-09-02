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

#include <RSAMPI/exanb_naive_grid.hxx>
#include <RSAMPI/exanb_naive_types.hxx>
#include <cstdint>
#include <median_of_medians.hxx>
#include <radius_generator.hxx>
#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbFinalizeCandidates : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT_OUTPUT, DocString{"Main grid that receives committed spheres"});
  ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{"Thickness (in cells) of the grid's ghost shell"});
  ADD_SLOT(GridT, candidate_grid, INPUT, REQUIRED,
           DocString{"Converged candidate grid (real cells = this rank's surviving candidates)"});
  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT_OUTPUT, DocString{"Radius generator state"});
  ADD_SLOT(uint64_t, next_id, INPUT_OUTPUT, uint64_t(0), DocString{"Next (rank-scoped) particle id to assign"});

  ADD_SLOT(uint64_t, nb_added_spheres, OUTPUT, DocString{"Number of spheres committed this round (global)"});
  ADD_SLOT(uint64_t, nb_added_spheres_local, OUTPUT, DocString{"Number of spheres committed this round (this rank only)"});

  inline void execute() override final {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const exanb::IJK dim = candidate_grid->dimension();
    const ssize_t gl = *ghost_layer;

    struct Survivor {
      exanb::Vec3d pos;
      double radius;
      int32_t phase;
      uint64_t priority;
    };
    std::vector<Survivor> survivors;
    for (ssize_t i = gl; i < dim.i - gl; i++) {
      for (ssize_t j = gl; j < dim.j - gl; j++) {
        for (ssize_t k = gl; k < dim.k - gl; k++) {
          auto& cell = candidate_grid->cell(exanb::IJK{i, j, k});
          const size_t n = cell.size();
          for (size_t s = 0; s < n; s++) {
            if (cell[exanb::field::confirmed][s] == 0) {
              continue;
            }
            if (cell[exanb::field::id][s] == exanb_naive::sentinel_id) {
              continue;
            }
            survivors.push_back(
                Survivor{exanb::Vec3d{cell[exanb::field::rx][s], cell[exanb::field::ry][s], cell[exanb::field::rz][s]},
                         cell[exanb::field::radius][s], cell[exanb::field::phase][s], cell[exanb::field::priority][s]});
          }
        }
      }
    }

    const uint64_t nb_spheres_total_max = RSARadiusGenerator->get_current_number();
    uint64_t local_count = survivors.size();
    uint64_t global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    if (global_count > nb_spheres_total_max) {
      std::vector<uint64_t> local_priorities(survivors.size());
      for (size_t i = 0; i < survivors.size(); i++) local_priorities[i] = survivors[i].priority;
      const double pivot = median_of_medians::find_pivot_for_first_elements(
          local_priorities, int64_t(nb_spheres_total_max), std::less<uint64_t>{});
      std::vector<Survivor> capped;
      capped.reserve(survivors.size());
      for (const auto& sv : survivors) {
        if (double(sv.priority) <= pivot) {
          capped.push_back(sv);
        }
      }
      survivors.swap(capped);
      local_count = survivors.size();
    }

    for (const auto& sv : survivors) {
      const uint64_t id = (uint64_t(rank) << 40) | (*next_id)++;
      exanb_naive::insert_sphere(*grid, sv.pos, id, sv.radius, sv.phase, sv.priority);
    }
    grid->rebuild_particle_offsets();

    uint64_t global_added = 0;
    MPI_Allreduce(&local_count, &global_added, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    *nb_added_spheres = global_added;
    *nb_added_spheres_local = local_count;
    RSARadiusGenerator->update_placed(global_added);
  }
};

template <class GridT>
using RSAMPIExanbFinalizeCandidatesTmpl = RSAMPIExanbFinalizeCandidates<GridT>;

ONIKA_AUTORUN_INIT(exanb_finalize_candidates) {
  OperatorNodeFactory::instance()->register_factory(
      "exanb_finalize_candidates", exanb::make_grid_variant_operator<RSAMPIExanbFinalizeCandidatesTmpl>);
}
}  // namespace rsa_mpi
