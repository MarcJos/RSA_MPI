#pragma once

#include <RSAMPI/fields.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>

#include <array>
#include <list_of_voxels.hxx>

namespace rsa_mpi {
namespace exanb_naive {

/// \brief Whether a voxel (origin_voxel + voxel_lengths, with corners relative
/// to origin_voxel in corners_voxel) is fully covered by a sphere in `grid`
/// enlarged by min_radius. IJK-native counterpart to
/// voxel_list::auxi::is_covered, scanning only the 3x3x3 neighborhood of the
/// voxel center's cell (real + ghost), like overlaps_existing.
template <class GridT>
inline bool voxel_is_covered(GridT& grid, const sac_de_billes::Point<3>& origin_voxel,
                             const sac_de_billes::Point<3>& voxel_lengths,
                             const std::array<sac_de_billes::Point<3>, 8>& corners_voxel, double min_radius) {
  const ::exanb::Vec3d center{origin_voxel[0] + 0.5 * voxel_lengths[0], origin_voxel[1] + 0.5 * voxel_lengths[1],
                              origin_voxel[2] + 0.5 * voxel_lengths[2]};
  const ::exanb::IJK loc = grid.locate_cell(center);
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
          const double sr = cell[::exanb::field::radius][s] + min_radius;
          const double sr2 = sr * sr;
          const double dx = cell[::exanb::field::rx][s] - origin_voxel[0];
          const double dy = cell[::exanb::field::ry][s] - origin_voxel[1];
          const double dz = cell[::exanb::field::rz][s] - origin_voxel[2];
          bool covers_all_corners = true;
          for (const auto& c : corners_voxel) {
            const double ddx = c[0] - dx;
            const double ddy = c[1] - dy;
            const double ddz = c[2] - dz;
            if (ddx * ddx + ddy * ddy + ddz * ddz > sr2) {
              covers_all_corners = false;
              break;
            }
          }
          if (covers_all_corners) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

/// \brief IJK-native counterpart to algorithm::auxi::update_covered_voxels:
/// once miss_rate exceeds desired_miss_rate, removes voxels fully covered by
/// `grid`'s spheres and subdivides the rest, using voxel_is_covered instead of
/// requiring an RSA_GRID. Leaves list_of_voxels itself untouched.
template <class GridT>
inline void update_covered_voxels(voxel_list::list_of_voxels<3>& uncovered_voxels, GridT& grid, double min_radius,
                                  double miss_rate, double desired_miss_rate) {
  if (miss_rate <= desired_miss_rate) {
    return;
  }
  auto is_covered = [&grid, min_radius](const sac_de_billes::Point<3>& origin_voxel,
                                        const sac_de_billes::Point<3>& voxel_lengths,
                                        const std::array<sac_de_billes::Point<3>, 8>& corners_voxel) {
    return voxel_is_covered(grid, origin_voxel, voxel_lengths, corners_voxel, min_radius);
  };
  uncovered_voxels.remove_covered_if(is_covered);
  uncovered_voxels.subdivide_uncovered_if(is_covered);
}

}  // namespace exanb_naive
}  // namespace rsa_mpi
