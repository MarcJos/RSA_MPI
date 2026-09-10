#pragma once

#include <exanb/amr/amr_grid.h>
#include <exanb/amr/amr_grid_algorithm.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <onika/log.h>

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace rsa_mpi {
namespace exanb_naive {

/// \brief Sorts `grid`'s particles into an adaptive per-cell sub-grid so
/// neighbor searches on crowded cells can window to nearby sub-cells instead
/// of scanning the whole cell. Reorders particles within each cell.
template <class GridT>
inline ::exanb::AmrGrid build_amr_index(GridT& grid, double avg_density = 2.0) {
  ::exanb::AmrGrid amr;
  ::exanb::rebuild_sub_grids(::onika::ldbg, grid, amr, avg_density);
  return amr;
}

/// \brief Calls f(cell, s) for particle slots in cell `nloc` within `reach`
/// of `query_p`. Falls back to a full scan if the cell has no sub-grid.
template <class GridT, class F>
inline void for_each_particle_near(GridT& grid, const ::exanb::AmrGrid& amr, const ::exanb::IJK& nloc, double reach,
                                   const ::exanb::Vec3d& query_p, F&& f) {
  auto& cell = grid.cell(nloc);
  const size_t n = cell.size();
  const size_t cell_i = grid.cell_index(nloc);
  const auto& sg_start = amr.sub_grid_start();
  if (sg_start.size() <= cell_i + 1) {
    for (size_t s = 0; s < n; s++) f(cell, s);
    return;
  }
  const size_t sgstart = sg_start[cell_i];
  const size_t sgsize = sg_start[cell_i + 1] - sgstart;
  if (sgsize == 0) {
    for (size_t s = 0; s < n; s++) f(cell, s);
    return;
  }

  const size_t n_sub_cells = sgsize + 1;
  const ssize_t sgside = static_cast<ssize_t>(std::floor(std::cbrt(double(n_sub_cells)) + 0.5));
  const ::exanb::Vec3d pc = grid.particle_pcoord(nloc, query_p);
  const double reach_frac = reach / grid.cell_size();

  const auto axis_range = [&](double c) -> std::pair<ssize_t, ssize_t> {
    ssize_t lo = static_cast<ssize_t>(std::floor((c - reach_frac) * sgside));
    ssize_t hi = static_cast<ssize_t>(std::floor((c + reach_frac) * sgside));
    lo = std::max<ssize_t>(lo, 0);
    hi = std::min<ssize_t>(hi, sgside - 1);
    return {lo, hi};
  };
  const auto [ilo, ihi] = axis_range(pc.x);
  const auto [jlo, jhi] = axis_range(pc.y);
  const auto [klo, khi] = axis_range(pc.z);
  if (ilo > ihi || jlo > jhi || klo > khi) {
    return;
  }

  const auto& sg_cells = amr.sub_grid_cells();
  for (ssize_t k = klo; k <= khi; k++) {
    for (ssize_t j = jlo; j <= jhi; j++) {
      for (ssize_t i = ilo; i <= ihi; i++) {
        const ssize_t sgindex = k * sgside * sgside + j * sgside + i;
        const size_t beginp = (sgindex == 0) ? 0 : sg_cells[sgstart + sgindex - 1];
        const size_t endp = (size_t(sgindex) < sgsize) ? sg_cells[sgstart + sgindex] : n;
        for (size_t s = beginp; s < endp; s++) f(cell, s);
      }
    }
  }
}

}  // namespace exanb_naive
}  // namespace rsa_mpi
