#include <mpi.h>
#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

// ExaNBody
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/make_grid_variant_operator.h>

// RSA MPI
#include <RSAMPI/fields.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace rsa_mpi {
using namespace onika;
using namespace onika::scg;

template <class GridT>
class RSAMPIExanbWriteXYZ : public OperatorNode {
  ADD_SLOT(GridT, grid, INPUT, REQUIRED, DocString{"exaNBody grid that contains particles"});
  ADD_SLOT(exanb::Domain, domain, INPUT, REQUIRED,
           DocString{"exaNBody domain (global bounds, written as the box size)"});
  ADD_SLOT(int, ghost_layer, INPUT, 1, DocString{"Thickness (in cells) of the grid's ghost shell"});
  ADD_SLOT(std::string, output_dir, INPUT, REQUIRED, DocString{"Root output directory"});
  ADD_SLOT(std::string, name, INPUT, std::string("packing.xyz"),
           DocString{"Output file name, written under output_dir"});

  inline void execute() override final {
    const exanb::IJK dim = grid->dimension();
    const ssize_t gl = *ghost_layer;

    std::ostringstream lines;
    size_t local_particles = 0;
    for (ssize_t i = gl; i < dim.i - gl; i++) {
      for (ssize_t j = gl; j < dim.j - gl; j++) {
        for (ssize_t k = gl; k < dim.k - gl; k++) {
          auto& cell = grid->cell(exanb::IJK{i, j, k});
          const size_t n = cell.size();
          for (size_t s = 0; s < n; s++) {
            lines << cell[exanb::field::type][s] << " " << cell[exanb::field::rx][s] << " "
                  << cell[exanb::field::ry][s] << " " << cell[exanb::field::rz][s] << "\n";
          }
          local_particles += n;
        }
      }
    }

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::string local_text = lines.str();
    const int local_bytes = static_cast<int>(local_text.size());
    const int local_count = static_cast<int>(local_particles);

    std::vector<int> byte_counts(size), particle_counts(size);
    MPI_Gather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_count, 1, MPI_INT, particle_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> displs(size, 0);
    int total_bytes = 0, total_particles = 0;
    if (rank == 0) {
      for (int p = 0; p < size; p++) {
        displs[p] = total_bytes;
        total_bytes += byte_counts[p];
        total_particles += particle_counts[p];
      }
    }

    std::vector<char> all_text(rank == 0 ? total_bytes : 0);
    MPI_Gatherv(local_text.data(), local_bytes, MPI_CHAR, all_text.data(), byte_counts.data(), displs.data(),
               MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      std::filesystem::create_directories(*output_dir);
      const std::string path = *output_dir + "/" + *name;
      std::ofstream out(path);
      if (!out) {
        onika::lerr << "Cannot open " << path << " for writing" << std::endl;
        return;
      }
      const exanb::Vec3d box = domain->bounds_size();
      out << total_particles << "\n" << box.x << " " << box.y << " " << box.z << "\n";
      std::stringstream body;
      body.write(all_text.data(), all_text.size());
      out << body.rdbuf();
    }
  }
};

template <class GridT>
using RSAMPIExanbWriteXYZTmpl = RSAMPIExanbWriteXYZ<GridT>;

ONIKA_AUTORUN_INIT(exanb_write_xyz) {
  OperatorNodeFactory::instance()->register_factory("exanb_write_xyz",
                                                    exanb::make_grid_variant_operator<RSAMPIExanbWriteXYZTmpl>);
}
}  // namespace rsa_mpi
