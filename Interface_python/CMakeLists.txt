add_subdirectory(pybind11)

find_package(OpenMP)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(pybind11 REQUIRED)

pybind11_add_module(rsa_mpi_py Interface.cxx)
target_link_libraries(rsa_mpi_py PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX interfMatools)
target_compile_definitions(rsa_mpi_py PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(rsa_mpi_py PUBLIC ../include/)

