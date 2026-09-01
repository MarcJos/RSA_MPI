#pragma once

#include <RSAMPI/fields.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace rsa_mpi {
namespace exanb_naive {

/// \brief Whether a candidate sphere (p, r) overlaps any particle already in `grid`.
/// Scans only the 3x3x3 neighborhood of p's cell (real + ghost), which is
/// sufficient since cell_size >= 2*max radius.
template <class GridT>
inline bool overlaps_existing(GridT& grid, const ::exanb::Vec3d& p, double r) {
  const ::exanb::IJK loc = grid.locate_cell(p);
  for (ssize_t di = -1; di <= 1; di++) {
    for (ssize_t dj = -1; dj <= 1; dj++) {
      for (ssize_t dk = -1; dk <= 1; dk++) {
        const ::exanb::IJK nloc{loc.i + di, loc.j + dj, loc.k + dk};
        if (!grid.contains(nloc)) {
          continue;
        }
        auto& cell = grid.cell(nloc);
        const size_t n = cell.size();
        for (size_t s = 0; s < n; s++) {
          const double dx = cell[::exanb::field::rx][s] - p.x;
          const double dy = cell[::exanb::field::ry][s] - p.y;
          const double dz = cell[::exanb::field::rz][s] - p.z;
          const double rr = cell[::exanb::field::radius][s] + r;
          if (dx * dx + dy * dy + dz * dz < rr * rr) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

/// \brief Like overlaps_existing, but only counts a neighbor with strictly
/// lower priority than `my_priority` (excluding `skip_id`).
/// \param require_confirmed if true, a neighbor only counts if its
/// `confirmed` field is also set (see resolve_candidates_pass).
template <class GridT>
inline bool overlaps_lower_priority(GridT& grid, const ::exanb::Vec3d& p, double r, uint64_t my_priority,
                                    uint64_t skip_id, bool require_confirmed = false) {
  const ::exanb::IJK loc = grid.locate_cell(p);
  for (ssize_t di = -1; di <= 1; di++) {
    for (ssize_t dj = -1; dj <= 1; dj++) {
      for (ssize_t dk = -1; dk <= 1; dk++) {
        const ::exanb::IJK nloc{loc.i + di, loc.j + dj, loc.k + dk};
        if (!grid.contains(nloc)) continue;
        auto& cell = grid.cell(nloc);
        const size_t n = cell.size();
        for (size_t s = 0; s < n; s++) {
          if (cell[::exanb::field::id][s] == skip_id) {
            continue;
          }
          if (cell[::exanb::field::priority][s] >= my_priority) {
            continue;
          }
          if (require_confirmed && cell[::exanb::field::confirmed][s] == 0) {
            continue;
          }
          const double dx = cell[::exanb::field::rx][s] - p.x;
          const double dy = cell[::exanb::field::ry][s] - p.y;
          const double dz = cell[::exanb::field::rz][s] - p.z;
          const double rr = cell[::exanb::field::radius][s] + r;
          if (dx * dx + dy * dy + dz * dz < rr * rr) return true;
        }
      }
    }
  }
  return false;
}

/// \brief One pass of the multi-rank candidate-resolution loop.
///
/// \return true if this rank changed anything this pass (caller
/// MPI_Allreduce-ORs this across ranks).
template <class GridT>
inline bool resolve_candidates_pass(GridT& main_grid, GridT& candidate_grid, ssize_t ghost_layer) {
  const ::exanb::IJK dim = candidate_grid.dimension();
  const ssize_t gl = ghost_layer;
  const size_t n_real_cells = size_t(dim.i - 2 * gl) * size_t(dim.j - 2 * gl) * size_t(dim.k - 2 * gl);

  // phase 1 (read-then-mutate, like the old single-phase version): decide,
  // from one consistent pre-pass snapshot, who dies to main_grid and who
  // becomes newly confirmed, then apply confirmations before phase 2 reads
  // them.
  std::vector<std::vector<size_t>> dead_per_cell;
  std::vector<std::vector<size_t>> newly_confirmed_per_cell;  // slot indices, within each cell
  dead_per_cell.reserve(n_real_cells);
  newly_confirmed_per_cell.reserve(n_real_cells);

  for (ssize_t i = gl; i < dim.i - gl; i++) {
    for (ssize_t j = gl; j < dim.j - gl; j++) {
      for (ssize_t k = gl; k < dim.k - gl; k++) {
        auto& cell = candidate_grid.cell(::exanb::IJK{i, j, k});
        const size_t n = cell.size();
        std::vector<size_t> dead;
        std::vector<size_t> newly_confirmed;
        for (size_t s = 0; s < n; s++) {
          if (cell[::exanb::field::confirmed][s] != 0) {
            continue;
          }
          const ::exanb::Vec3d p{cell[::exanb::field::rx][s], cell[::exanb::field::ry][s], cell[::exanb::field::rz][s]};
          const double r = cell[::exanb::field::radius][s];
          const uint64_t id = cell[::exanb::field::id][s];
          const uint64_t priority = cell[::exanb::field::priority][s];
          if (overlaps_existing(main_grid, p, r)) {
            dead.push_back(s);
          } else if (!overlaps_lower_priority(candidate_grid, p, r, priority, id)) {
            newly_confirmed.push_back(s);
          }
        }
        dead_per_cell.push_back(std::move(dead));
        newly_confirmed_per_cell.push_back(std::move(newly_confirmed));
      }
    }
  }

  bool changed = false;
  {
    size_t cell_idx = 0;
    for (ssize_t i = gl; i < dim.i - gl; i++)
      for (ssize_t j = gl; j < dim.j - gl; j++)
        for (ssize_t k = gl; k < dim.k - gl; k++) {
          auto& cell = candidate_grid.cell(::exanb::IJK{i, j, k});
          for (size_t s : newly_confirmed_per_cell[cell_idx]) {
            cell[::exanb::field::confirmed][s] = 1;
            changed = true;
          }
          cell_idx++;
        }
  }

  // phase 2: now that this pass's confirmations are visible everywhere
  // (including to neighbors in earlier-visited cells), reject any remaining
  // not-yet-confirmed survivor that overlaps a CONFIRMED lower-priority
  // neighbor.
  {
    std::vector<bool> already_dead;
    size_t cell_idx = 0;
    for (ssize_t i = gl; i < dim.i - gl; i++)
      for (ssize_t j = gl; j < dim.j - gl; j++)
        for (ssize_t k = gl; k < dim.k - gl; k++) {
          auto& cell = candidate_grid.cell(::exanb::IJK{i, j, k});
          auto& dead = dead_per_cell[cell_idx++];
          const size_t n = cell.size();
          already_dead.clear();
          already_dead.resize(n);  // values are set to false
          for (size_t s : dead) {
            already_dead[s] = true;
          }
          for (size_t s = 0; s < n; s++) {  // already dead from phase 1 (main_grid).
            if (already_dead[s]) {
              continue;
            }
            if (cell[::exanb::field::confirmed][s] != 0) {  // confirmed this pass or earlier - safe.
              continue;
            }
            const ::exanb::Vec3d p{cell[::exanb::field::rx][s], cell[::exanb::field::ry][s],
                                   cell[::exanb::field::rz][s]};
            const double r = cell[::exanb::field::radius][s];
            const uint64_t id = cell[::exanb::field::id][s];
            const uint64_t priority = cell[::exanb::field::priority][s];
            if (overlaps_lower_priority(candidate_grid, p, r, priority, id, true)) {
              dead.push_back(s);
            }
          }
        }
  }

  size_t cell_idx = 0;
  for (ssize_t i = gl; i < dim.i - gl; i++) {
    for (ssize_t j = gl; j < dim.j - gl; j++) {
      for (ssize_t k = gl; k < dim.k - gl; k++) {
        auto& dead = dead_per_cell[cell_idx++];
        if (dead.empty()) continue;
        std::sort(dead.begin(), dead.end());
        dead.erase(std::unique(dead.begin(), dead.end()), dead.end());
        changed = true;
        auto& cell = candidate_grid.cell(::exanb::IJK{i, j, k});
        // remove highest slot index first so earlier removals in this cell
        // aren't invalidated by the swap-with-last-and-shrink below.
        for (auto it = dead.rbegin(); it != dead.rend(); ++it) {
          const size_t s = *it;
          const size_t last = cell.size() - 1;
          if (s != last) cell.swap(s, last);
          cell.resize(last, candidate_grid.cell_allocator());
        }
      }
    }
  }
  return changed;
}

}  // namespace exanb_naive
}  // namespace rsa_mpi
