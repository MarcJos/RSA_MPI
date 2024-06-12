//! Copyright : Apache 2.0, see LICENSE 
//! 
#include<test.hxx>
#include<chrono>
#include <omp.h>
#include <rsa_parameters.hxx>
#include<rsa_test_perf.hxx>
//#define VERBOSE_FILL 1

int main(int argc, char** argv) {
  rsa_mpi::init(&argc, &argv);

  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int nb_threads = 1;

  rsa_mpi::rsa_parameters param = rsa_mpi::read_input(argc, argv);
  param.minimal_requirement();
  param.display();

  auto test_gen = [&](rsa_mpi::rsa_parameters& p) {
    test_perf(p);
    return true;
    };

  make_XP(nb_threads, test_gen, param);

  rsa_mpi::finalize();
  return EXIT_SUCCESS;
}
