#include <mpi.h>
#include <onika/log.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

// RSA MPI
#include <RSAMPI/fields.h>

#include <RSAMPI/exanb_naive_draw.hxx>
#include <RSAMPI/exanb_naive_types.hxx>
#include <radius_generator.hxx>
#include <random>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbNaiveCandidates : public OperatorNode {
  ADD_SLOT(exanb::Vec3d, inf, INPUT, REQUIRED, DocString{"Minimum coordinates of the system"});
  ADD_SLOT(exanb::Vec3d, sup, INPUT, REQUIRED, DocString{"Maximum coordinates of the system"});
  ADD_SLOT(int, size, INPUT, REQUIRED, DocString{"Number of candidates drawn per iteration"});
  ADD_SLOT(int, seed, INPUT, 0, DocString{"RNG seed"});

  ADD_SLOT(sac_de_billes::RadiusGenerator<3>, RSARadiusGenerator, INPUT_OUTPUT, DocString{"Radius generator state"});
  ADD_SLOT(std::mt19937, rng, INPUT_OUTPUT, std::mt19937(), DocString{"RNG state, kept across iterations"});
  ADD_SLOT(bool, rng_initialized, INPUT_OUTPUT, false, DocString{"Whether rng has been seeded yet"});

  ADD_SLOT(exanb_naive::Candidates, candidates, OUTPUT, DocString{"Drawn candidates for this iteration"});

  inline void execute() override final {
    if (!*rng_initialized) {
      rng->seed(*seed);
      *rng_initialized = true;
    }
    *candidates = exanb_naive::draw_candidates(*inf, *sup, *size, *RSARadiusGenerator, *rng);
  }
};

template <class GridT>
using RSAMPIExanbNaiveCandidatesTmpl = RSAMPIExanbNaiveCandidates<GridT>;

ONIKA_AUTORUN_INIT(exanb_naive_candidates) {
  OperatorNodeFactory::instance()->register_factory("exanb_naive_candidates",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbNaiveCandidatesTmpl>);
}
}  // namespace rsa_mpi
