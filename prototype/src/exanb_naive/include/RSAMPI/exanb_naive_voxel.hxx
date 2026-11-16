#pragma once

#include <RSAMPI/exanb_naive_amr.hxx>
#include <RSAMPI/fields.h>
#include <exanb/amr/amr_grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>

#include <array>
#include <cmath>
#include <list_of_voxels.hxx>

namespace rsa_mpi {
namespace exanb_naive {

/// \brief Whether a voxel (origin_voxel + voxel_lengths, with corners relative
/// to origin_voxel in corners_voxel) is fully covered by a sphere in `grid`
/// enlarged by min_radius. IJK-native counterpart to
/// voxel_list::auxi::is_covered, scanning only the 3x3x3 neighborhood of the
/// voxel center's cell (real + ghost), like overlaps_existing; within each
/// neighbor cell, `amr` narrows the scan to the sub-cells actually within reach.
template <class GridT>
inline bool voxel_is_covered(GridT& grid, const ::exanb::AmrGrid& amr, const sac_de_billes::Point<3>& origin_voxel,
                             const sac_de_billes::Point<3>& voxel_lengths,
                             const std::array<sac_de_billes::Point<3>, 8>& corners_voxel, double min_radius,
                             double max_radius) {
  const ::exanb::Vec3d center{origin_voxel[0] + 0.5 * voxel_lengths[0], origin_voxel[1] + 0.5 * voxel_lengths[1],
                              origin_voxel[2] + 0.5 * voxel_lengths[2]};
  const ::exanb::IJK loc = grid.locate_cell(center);
  const double half_diagonal =
      0.5 * std::sqrt(voxel_lengths[0] * voxel_lengths[0] + voxel_lengths[1] * voxel_lengths[1] +
                      voxel_lengths[2] * voxel_lengths[2]);
  const double reach = half_diagonal + max_radius + min_radius;
  bool found = false;
  for (ssize_t di = -1; di <= 1 && !found; di++) {
    for (ssize_t dj = -1; dj <= 1 && !found; dj++) {
      for (ssize_t dk = -1; dk <= 1 && !found; dk++) {
        const ::exanb::IJK nloc{loc.i + di, loc.j + dj, loc.k + dk};
        if (!grid.contains(nloc)) {
          continue;
        }
        for_each_particle_near(grid, amr, nloc, reach, center, [&](auto& cell, size_t s) {
          if (found) return;
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
            found = true;
          }
        });
      }
    }
  }
  return found;
}

/// \brief Local (this rank's) outcome of a call to update_covered_voxels.
struct VoxelRefinementStats {
  bool refined = false;   ///< true if remove_covered_if/subdivide_uncovered_if ran this call
  int64_t covered = 0;    ///< voxels found fully covered and removed
  int64_t generated = 0;  ///< voxel count after subdivision
};

/// \brief IJK-native counterpart to algorithm::auxi::update_covered_voxels:
/// once miss_rate exceeds desired_miss_rate, removes voxels fully covered by
/// `grid`'s spheres and subdivides the rest, using voxel_is_covered instead of
/// requiring an RSA_GRID. Leaves list_of_voxels itself untouched.
template <class GridT>
inline VoxelRefinementStats update_covered_voxels(voxel_list::list_of_voxels<3>& uncovered_voxels, GridT& grid,
                                                  double min_radius, double max_radius, double miss_rate,
                                                  double desired_miss_rate) {
  if (miss_rate <= desired_miss_rate) {
    return VoxelRefinementStats{};
  }
  const ::exanb::AmrGrid amr = build_amr_index(grid);
  auto is_covered = [&grid, &amr, min_radius, max_radius](const sac_de_billes::Point<3>& origin_voxel,
                                                          const sac_de_billes::Point<3>& voxel_lengths,
                                                          const std::array<sac_de_billes::Point<3>, 8>& corners_voxel) {
    return voxel_is_covered(grid, amr, origin_voxel, voxel_lengths, corners_voxel, min_radius, max_radius);
  };
  const int64_t before = int64_t(uncovered_voxels.size());
  uncovered_voxels.remove_covered_if(is_covered);
  const int64_t after_remove = int64_t(uncovered_voxels.size());
  uncovered_voxels.subdivide_uncovered_if(is_covered);
  return VoxelRefinementStats{true, before - after_remove, int64_t(uncovered_voxels.size())};
}

}  // namespace exanb_naive
}  // namespace rsa_mpi
