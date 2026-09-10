set(rsa_mpi_exanb_naive_INCLUDE_DIRS ${RSA_MPI_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(rsa_mpi_exanb_naive_LINK_LIBRARIES exanbCore exanbGridCellParticles exanbAMR)
onika_add_plugin(rsa_mpi_exanb_naive ${CMAKE_CURRENT_SOURCE_DIR})
