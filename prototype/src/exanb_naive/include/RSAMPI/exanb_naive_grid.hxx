#pragma once

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>

#include <RSAMPI/exanb_naive_types.hxx>
#include <cstdint>

namespace rsa_mpi {
namespace exanb_naive {

/// \brief Inserts a sphere into `grid`'s owning cell (no ghost propagation).
template <class GridT>
inline void insert_sphere(GridT& grid, const ::exanb::Vec3d& p, uint64_t id, double r, int32_t phase, uint64_t priority,
                          int32_t confirmed = 0) {
  const ::exanb::IJK loc = grid.locate_cell(p);
  ParticleTuple pt(p.x, p.y, p.z, id, r, phase, priority, confirmed);
  grid.cell(loc).push_back(pt, grid.cell_allocator());
}

/// \brief Builds an empty grid with the same geometry as `grid` (including
/// max_neighbor_distance, needed for ghost_comm_scheme to compute the
/// correct ghost_layers()).
template <class GridT>
inline GridT make_scratch_grid(const GridT& grid) {
  GridT scratch;
  scratch.set_origin(grid.origin());
  scratch.set_cell_size(grid.cell_size());
  scratch.set_offset(grid.offset());
  scratch.set_dimension(grid.dimension());
  scratch.set_max_neighbor_distance(grid.max_neighbor_distance());
  return scratch;
}

}  // namespace exanb_naive
}  // namespace rsa_mpi
