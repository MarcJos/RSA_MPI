#pragma once

#include <RSAMPI/fields.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid_fields.h>
#include <onika/soatl/field_tuple.h>

#include <cstdint>
#include <limits>
#include <vector>

namespace rsa_mpi {
namespace exanb_naive {

/// One round's drawn candidate spheres (position, radius, type, priority),
/// in struct-of-arrays form.
struct Candidates {
  std::vector<::exanb::Vec3d> pos;
  std::vector<double> radius;
  std::vector<::exanb::ParticleTypeInt> type;
  std::vector<uint64_t> priority;
  size_t size() const { return pos.size(); }
};

/// Per-particle field tuple stored in a candidate_grid/main grid cell.
using ParticleTuple = ::onika::soatl::FieldTuple<::exanb::field::_rx, ::exanb::field::_ry, ::exanb::field::_rz,
                                                 ::exanb::field::_id, ::exanb::field::_radius, ::exanb::field::_type,
                                                 ::exanb::field::_priority, ::exanb::field::_confirmed>;

inline constexpr uint64_t sentinel_id = std::numeric_limits<uint64_t>::max();

}  // namespace exanb_naive
}  // namespace rsa_mpi
