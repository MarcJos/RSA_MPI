//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <stdlib.h>

namespace rsa_mpi {
/*
 * struct rsa_parameters
 *
 * A struct to hold parameters for RSA-based parallel processing.
 */
struct rsa_parameters {
  int DIM = -1;                               // Dimension of the system
  double radius = -1;                         // Radius of particles
  std::vector<double> l_min = { 0,0,0 };        // Minimum coordinates of the system
  std::vector<double> l_max = { 1,1,1 };        // Maximum coordinates of the system
  int seed = 0;                               // Seed for random number generation
  int n_draw = 10;                            // Number of draws
  int size = 6000;                            // TO DO RENAME
  bool paraview = 0;                           // Flag for Paraview visualization


  /*
   * Displays the values of RSA parameters.
   */
  void display() {
    int mpi_rank = rsa_mpi::get_my_rank();
    if (mpi_rank == 0) {
      std::cout << "============= Parameters =================" << std::endl;
      std::cout << "Dim         = " << DIM << std::endl;
      std::cout << "Radius      = " << radius << std::endl;
      std::cout << "Domain size = ";
      for (int d = 0; d < DIM - 1; d++) std::cout << "[" << l_min[d] << "," << l_max[d] << "]x";
      std::cout << "[" << l_min[DIM - 1] << "," << l_max[DIM - 1] << "]" << std::endl;
      std::cout << "Seed        = " << seed << std::endl;
      std::cout << "N Draw      = " << n_draw << std::endl;
      std::cout << "Size        = " << size << std::endl;
      std::cout << "Paraview    = " << paraview << " (0 disabled, 1 activated)" << std::endl;
      std::cout << "==========================================" << std::endl;
    }
  }

  void minimal_requirement() {
    bool okay = true;
    int mpi_rank = rsa_mpi::get_my_rank();
    if (mpi_rank == 0) {
      if (DIM < 2) {
        std::cout << "DIM is not defined correctly, Value=" << DIM << "." << std::endl;
        std::cout << "Please use --dim Value." << std::endl;
        okay = false;
      }
      if (radius < 0.0) {
        std::cout << "Radius is not defined correctly, Value=" << radius << "." << std::endl;
        std::cout << "Please use --radius Value." << std::endl;
        okay = false;
      }

      if (l_min.size() != DIM) { okay = false; std::cout << "Domain inf is not defined correctly, DIM= " << l_min.size() << " instead of " << DIM << "." << std::endl; }
      if (l_max.size() != DIM) { okay = false; std::cout << "Domain inf is not defined correctly, DIM= " << l_min.size() << " instead of " << DIM << "." << std::endl; }

      for (int d = 0; d < DIM; d++) {
        if (l_min[d] >= l_max[d]) {
          okay = false; std::cout << "The domain simulation is not correctly defined, for DIM= " << d << " , inf >= max." << std::endl;
        }
      }
      if (!okay) std::abort;
    }

  };
};

/*
 * Display help command.
 *
*/
inline void help()
{
    if (rsa_mpi::get_my_rank() == 0) {
      std::cout << "=== Available options ===" << std::endl;
      std::cout << "--dim [int, >= 2],           set dimension." << std::endl;
      std::cout << "--radius [double],           set sphere radius." << std::endl;
      std::cout << "--seed [int],                set seed used to generate random numbers." << std::endl;
      std::cout << "--size [int],                set the maximum number of spheres thrown per dart." << std::endl;
      std::cout << "--inf [double, ..., double], set the position of the bottom limit of the simulation domain." << std::endl;
      std::cout << "--sup [double, ..., double], set the position of the upper limit of the simulation domain." << std::endl;
      std::cout << "--paraview [int, 0 or 1],    disable [0] or activate [1] the feature: write paraview files, default is 0." << std::endl;
    }
}

/*
 * Reads RSA parameters from command line arguments.
 *
 * argc: Number of command line arguments.
 * argv: Array of command line argument strings.
 *
 * Returns: An rsa_parameters struct with the parsed parameters.
 */
rsa_parameters read_input(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Wrong number of arguments" << std::endl;
    help();
    std::exit(EXIT_FAILURE);
  }

  if (std::string(argv[1]) == "--help") {
    help();
    std::exit(EXIT_FAILURE);
  }


  rsa_parameters param;
  int i = 1;
  int dim_l = 0;

  while (i < argc) {
    std::string key = argv[i];
    if (i == 1 && key != "--dim") {
      std::cout << "Wrong arg: " << argv[i] << std::endl;
      std::cout << "Please, use the dimension input parameter as first input parameter (--dim)" << std::endl;
      std::abort();
    } else if (i == 1) {
      param.DIM = atoi(argv[i + 1]);
      i += 2;
      param.l_min.assign(param.DIM, 0);
      param.l_max.assign(param.DIM, 0);
      continue;
    }

    if (key == "--radius") {
      param.radius = atof(argv[i + 1]);
      i += 2;
      continue;
    }

    if (key == "--seed") {
      param.seed = atoi(argv[i + 1]);
      i += 2;
      continue;
    }

    if (key == "--n_draw") {
      param.n_draw = atoi(argv[i + 1]);
      i += 2;
      continue;
    }

    if (key == "--size") {
      param.size = atoi(argv[i + 1]);
      i += 2;
      continue;
    }

    if (key == "--inf") {
      for (int d = 0; d < param.DIM; d++) {
        param.l_min[d] = atof(argv[i + 1 + d]);
      }
      i += 1 + param.DIM;
      continue;
    }

    if (key == "--sup") {
      for (int d = 0; d < param.DIM; d++) {
        param.l_max[d] = atof(argv[i + 1 + d]);
      }
      i += 1 + param.DIM;
      continue;
    }

    if (key == "--paraview") {
      param.paraview = atoi(argv[i + 1]);
      i += 2;
      continue;
    }
    // no key found
    if (rsa_mpi::get_my_rank() == 0)  std::cout << key << " is ignored." << std::endl;
    i++;
  }

  if (i != argc) {
    std::cout << "Wrong number of arguments" << std::endl;
    std::abort();
  }

  param.minimal_requirement();
  return param;
}
}
