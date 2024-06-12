add_library(rsa_mpi_lib operator_algorithm.cxx)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_include_directories(rsa_mpi_lib PUBLIC ../include/)
target_link_libraries(rsa_mpi_lib MPI::MPI_CXX interfMatools)
