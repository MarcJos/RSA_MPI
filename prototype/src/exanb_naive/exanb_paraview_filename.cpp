#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/string_utils.h>

#include <string>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

class RSAMPIExanbParaviewFilename : public OperatorNode {
  ADD_SLOT(std::string, output_dir, INPUT, REQUIRED, DocString{"Root output directory"});
  ADD_SLOT(std::string, name, INPUT, std::string("rsa_mpi"), DocString{"Base name of the Paraview output files"});
  ADD_SLOT(int, iteration, PRIVATE, 0);
  ADD_SLOT(std::string, filename, OUTPUT, DocString{"Path prefix (no extension) for this round's write_paraview dump"});

  inline void execute() override final {
    ++(*iteration);
    *filename = onika::format_string("%s/ParaviewOutputDir/%s_%05d", *output_dir, *name, *iteration);
  }
};

class RSAMPIExanbParaviewFinalFilename : public OperatorNode {
  ADD_SLOT(std::string, output_dir, INPUT, REQUIRED, DocString{"Root output directory"});
  ADD_SLOT(std::string, name, INPUT, std::string("rsa_mpi"), DocString{"Base name of the Paraview output files"});
  ADD_SLOT(std::string, filename, OUTPUT, DocString{"Path prefix (no extension) for the final write_paraview dump"});

  inline void execute() override final { *filename = *output_dir + "/ParaviewOutputDir/" + *name + "_final"; }
};

ONIKA_AUTORUN_INIT(exanb_paraview_filename) {
  OperatorNodeFactory::instance()->register_factory("exanb_paraview_filename",
                                                    make_simple_operator<RSAMPIExanbParaviewFilename>);
  OperatorNodeFactory::instance()->register_factory("exanb_paraview_final_filename",
                                                    make_simple_operator<RSAMPIExanbParaviewFinalFilename>);
}
}  // namespace rsa_mpi
