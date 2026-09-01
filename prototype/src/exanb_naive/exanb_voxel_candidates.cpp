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
#include <RSAMPI/exanb_naive_draw.hxx>
#include <RSAMPI/exanb_naive_types.hxx>
#include <list_of_voxels.hxx>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>
#include <random>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbVoxelCandidates : public OperatorNode {
  ADD_SLOT(voxel_list::list_of_voxels<3>, uncovered_voxels, INPUT, REQUIRED,
           DocString{"Voxels not yet covered by any placed sphere"});
  ADD_SLOT(int, size, INPUT, REQUIRED,
           DocString{"Ceiling on the (expected) number of shots per draw, across all MPI processes"});
  ADD_SLOT(int, seed, INPUT, 0, DocString{"RNG seed"});

  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT_OUTPUT, DocString{"Radius generator state"});
  ADD_SLOT(std::mt19937, rng, INPUT_OUTPUT, std::mt19937(), DocString{"RNG state, kept across iterations"});
  ADD_SLOT(bool, rng_initialized, INPUT_OUTPUT, false, DocString{"Whether rng has been seeded yet"});

  ADD_SLOT(exanb_naive::Candidates, candidates, OUTPUT, DocString{"Drawn candidates for this iteration"});
  ADD_SLOT(int, shots, OUTPUT, DocString{"Number of shots drawn (locally, Poisson-distributed) this iteration"});

  inline void execute() override final {
    if (!*rng_initialized) {
      rng->seed(*seed);
      *rng_initialized = true;
    }

    const double intensity = algorithm::auxi::compute_intensity_poisson(int64_t(uncovered_voxels->size()));
    double maximum_intensity = 0.0;
    MPI_Allreduce(&intensity, &maximum_intensity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    *shots = algorithm::auxi::compute_nb_shots_voxel(*rng, intensity, maximum_intensity, *size);

    const auto positions = uncovered_voxels->pick_points(*shots, *rng);
    *candidates = exanb_naive::draw_candidates_at(positions, *RSARadiusGenerator, *rng);
  }
};

template <class GridT>
using RSAMPIExanbVoxelCandidatesTmpl = RSAMPIExanbVoxelCandidates<GridT>;

ONIKA_AUTORUN_INIT(exanb_voxel_candidates) {
  OperatorNodeFactory::instance()->register_factory("exanb_voxel_candidates",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbVoxelCandidatesTmpl>);
}
}  // namespace rsa_mpi
