#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <iomanip>
#include <sstream>
#include <string>

// RSA MPI
#include <rsa_domain.hxx>
#include <rsa_grid.hxx>

//
#include <operator_paraview.hxx>

// RSA MPI
#include <RSAMPI/make_variant_operator.hpp>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

// subdirectory created inside output_dir, holding every Paraview output of a run
static const char* PARAVIEW_SUBDIR = "ParaviewOutputDir";

template <int DIM>
class RSAMPIParaview : public OperatorNode {
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, REQUIRED, DocString{"Grid that contains particles"});
  ADD_SLOT(std::string, name, INPUT, std::string("rsa_mpi"), DocString{"Base name of the Paraview output files"});
  ADD_SLOT(std::string, output_dir, INPUT, REQUIRED,
           DocString{"Root output directory; Paraview files are written under <output_dir>/ParaviewOutputDir/. "
                     "Defined by default in rsa_mpi_default.msp (global: output_dir)."});
  ADD_SLOT(int, iteration, PRIVATE, 0);

  inline std::string documentation() const override final {
    return R"EOF(
Writes the current grid content to a Paraview (.pvtp/.vtp) file, one call per
algorithm iteration. Files are numbered "<name>_00001.pvtp", "<name>_00002.pvtp", ...
and written under "<output_dir>/ParaviewOutputDir/".

      YAML example:

        naive_algorithm:
          loop: true
          condition: continue_loop
          body:
            - naive_draw
            - rsa_mpi_paraview
            - naive_draw_continue
      )EOF";
  }

  inline void execute() override final {
    ++(*iteration);
    std::ostringstream oss;
    oss << *name << "_" << std::setw(5) << std::setfill('0') << *iteration;

    const std::string directory = *output_dir + "/" + PARAVIEW_SUBDIR;
    rsa_paraview::paraview<DIM>(directory, oss.str(), *RSAGrid);
  }
};

template <int DIM>
class RSAMPIParaviewFinal : public OperatorNode {
  ADD_SLOT(rsa_grid<DIM>, RSAGrid, INPUT, REQUIRED, DocString{"Grid that contains particles"});
  ADD_SLOT(std::string, name, INPUT, std::string("rsa_mpi"), DocString{"Base name of the Paraview output files"});
  ADD_SLOT(std::string, output_dir, INPUT, REQUIRED,
           DocString{"Root output directory; Paraview files are written under <output_dir>/ParaviewOutputDir/. "
                     "Defined by default in rsa_mpi_default.msp (global: output_dir)."});

  inline std::string documentation() const override final {
    return R"EOF(
Writes the final grid content to a Paraview (.pvtp/.vtp) file, named
"<name>_final.pvtp", under "<output_dir>/ParaviewOutputDir/". Meant to be
called once, at the end of the simulation.

      YAML example:

        simulation_epilog:
          - rsa_mpi_paraview_final
      )EOF";
  }

  inline void execute() override final {
    const std::string directory = *output_dir + "/" + PARAVIEW_SUBDIR;
    const std::string basename = *name + "_final";
    rsa_paraview::paraview<DIM>(directory, basename, *RSAGrid);
  }
};

template <int DIM>
using RSAMPIParaviewTmpl = RSAMPIParaview<DIM>;
template <int DIM>
using RSAMPIParaviewFinalTmpl = RSAMPIParaviewFinal<DIM>;

ONIKA_AUTORUN_INIT(paraview) {
  OperatorNodeFactory::instance()->register_factory("rsa_mpi_paraview", make_rsa_mpi_operator<RSAMPIParaviewTmpl>);
  OperatorNodeFactory::instance()->register_factory("rsa_mpi_paraview_final",
                                                    make_rsa_mpi_operator<RSAMPIParaviewFinalTmpl>);
}
}  // namespace rsa_mpi
