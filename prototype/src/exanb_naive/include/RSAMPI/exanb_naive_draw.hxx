#pragma once

#include <exanb/core/domain.h>

#include <RSAMPI/exanb_naive_types.hxx>
#include <basic_types.hxx>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>
#include <random>
#include <vector>

namespace rsa_mpi {
namespace exanb_naive {

/// \brief Draws `n` candidate spheres uniformly in [inf,sup].
/// \param radius_gen source of per-candidate radius/type.
/// \param rng consumed for position, priority, radius and type, in that order.
inline Candidates draw_candidates(const ::exanb::Vec3d& inf, const ::exanb::Vec3d& sup, int n,
                                  sac_de_billes::RadiusGenerator<3>& radius_gen, std::mt19937& rng) {
  Candidates c;
  c.pos.resize(n);
  {
    std::uniform_real_distribution<double> law(inf.x, sup.x);
    for (int i = 0; i < n; i++) {
      c.pos[i].x = law(rng);
    }
  }
  {
    std::uniform_real_distribution<double> law(inf.y, sup.y);
    for (int i = 0; i < n; i++) {
      c.pos[i].y = law(rng);
    }
  }
  {
    std::uniform_real_distribution<double> law(inf.z, sup.z);
    for (int i = 0; i < n; i++) {
      c.pos[i].z = law(rng);
    }
  }

  auto priorities = algorithm::generate_priority<int>(n, rng);
  c.priority.assign(priorities.begin(), priorities.end());

  auto phases_radii = radius_gen(n, rng);
  const auto& phases = std::get<0>(phases_radii);
  const auto& radii = std::get<1>(phases_radii);
  c.radius.assign(radii.begin(), radii.end());
  c.type.assign(phases.begin(), phases.end());

  return c;
}

/// \brief Like draw_candidates, but positions are given (e.g. picked from a
/// voxel_list::list_of_voxels) instead of drawn uniformly in [inf,sup].
inline Candidates draw_candidates_at(const std::vector<sac_de_billes::Point<3>>& positions,
                                     sac_de_billes::RadiusGenerator<3>& radius_gen, std::mt19937& rng) {
  Candidates c;
  const int n = static_cast<int>(positions.size());
  c.pos.resize(n);
  for (int i = 0; i < n; i++) {
    c.pos[i] = ::exanb::Vec3d{positions[i][0], positions[i][1], positions[i][2]};
  }

  auto priorities = algorithm::generate_priority<int>(n, rng);
  c.priority.assign(priorities.begin(), priorities.end());

  auto phases_radii = radius_gen(n, rng);
  const auto& phases = std::get<0>(phases_radii);
  const auto& radii = std::get<1>(phases_radii);
  c.radius.assign(radii.begin(), radii.end());
  c.type.assign(phases.begin(), phases.end());

  return c;
}

}  // namespace exanb_naive
}  // namespace rsa_mpi
