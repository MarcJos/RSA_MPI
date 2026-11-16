#include <mpi.h>
#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/grid_cell_particles/particle_region.h>

// RSA MPI
#include <RSAMPI/fields.h>

#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbRegionFilter : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT_OUTPUT, DocString{"Main grid to filter in place (real cells only)"});
  ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{"Thickness (in cells) of the grid's ghost shell"});
  ADD_SLOT(exanb::ParticleRegions, particle_regions, INPUT, OPTIONAL,
           DocString{"Named regions (bounds/quadric), referenced by name in region"});
  ADD_SLOT(exanb::ParticleRegionCSG, region, INPUT, OPTIONAL,
           DocString{"Boolean expression of region names; particles outside it are removed. No-op if unset."});

  inline void execute() override final {
    if (!region.has_value()) {
      return;
    }
    if (!particle_regions.has_value()) {
      fatal_error() << "region is defined, but particle_regions has no value" << std::endl;
    }
    if (region->m_nb_operands == 0) {
      region->build_from_expression_string(particle_regions->data(), particle_regions->size());
    }
    const exanb::ParticleRegionCSGShallowCopy prcsg = *region;

    const exanb::IJK dim = grid->dimension();
    const ssize_t gl = *ghost_layer;

    uint64_t local_removed = 0;
    std::vector<size_t> dead;
    for (ssize_t i = gl; i < dim.i - gl; i++) {
      for (ssize_t j = gl; j < dim.j - gl; j++) {
        for (ssize_t k = gl; k < dim.k - gl; k++) {
          auto& cell = grid->cell(exanb::IJK{i, j, k});
          const size_t n = cell.size();
          dead.clear();
          for (size_t s = 0; s < n; s++) {
            const exanb::Vec3d p{cell[exanb::field::rx][s], cell[exanb::field::ry][s], cell[exanb::field::rz][s]};
            if (!prcsg.contains(p, cell[exanb::field::id][s])) {
              dead.push_back(s);
            }
          }
          if (dead.empty()) {
            continue;
          }
          local_removed += dead.size();
          for (auto it = dead.rbegin(); it != dead.rend(); ++it) {
            const size_t s = *it;
            const size_t last = cell.size() - 1;
            if (s != last) cell.swap(s, last);
            cell.resize(last, grid->cell_allocator());
          }
        }
      }
    }
    grid->rebuild_particle_offsets();

    uint64_t global_removed = 0;
    MPI_Allreduce(&local_removed, &global_removed, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      onika::lout << "exanb_region_filter: removed " << global_removed << " particle(s) outside region" << std::endl;
    }
  }
};

template <class GridT>
using RSAMPIExanbRegionFilterTmpl = RSAMPIExanbRegionFilter<GridT>;

ONIKA_AUTORUN_INIT(exanb_region_filter) {
  OperatorNodeFactory::instance()->register_factory("exanb_region_filter",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbRegionFilterTmpl>);
}
}  // namespace rsa_mpi
