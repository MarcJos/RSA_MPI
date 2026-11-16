
project(
  onika-mpi-rsa
  VERSION 0.0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===========================
# === CMake customization ===
# ===========================
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release RelWithDebInfo Debug)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

## Onika
find_package(OpenMP REQUIRED)
find_package(onika REQUIRED)

## exaNBody
# these must be set before find_package(exaNBody): they configure how exaNBody
# builds its grid-related targets (domain-specific particle fields, in particular)
set(XNB_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XNB_COMPILE_DEFINITIONS "XNB_DOMAIN_SPECIFIC_FIELDS_INCLUDE=<RSAMPI/fields.h>" RSA_MPI_VERSION=${PROJECT_VERSION})
set(XNB_APP_HAS_GRID_SPECIALIZATIONS ON)

# Only exanbGridCellParticles (init_rcb_grid), exanbMPI (ghost_comm_scheme,
# exanb_update_ghosts) and exanbIO (write_paraview) are actually used by our
# operators; the rest are MD/DEM-specific plugins we never reference.
foreach(_xnb_unused_plugin exanbAMR exanbAnalytics exanbCompute exanbDebug exanbDefBox
                           exanbExtraStorage exanbLogic exanbMesh exanbParticleNeighbors)
  set(ONIKA_BUILD_${_xnb_unused_plugin} OFF CACHE BOOL "Build ${_xnb_unused_plugin}" FORCE)
endforeach()

find_package(exaNBody REQUIRED)

set(rsa_mpi_DIR ${CMAKE_SOURCE_DIR})
message(STATUS "${rsa_mpi_DIR}/include")

# ======================================
# === add compiler option with CUDA  ===
# ======================================
if(ONIKA_BUILD_CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
endif()


if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS_DEBUG   "-g -traceback")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHOST")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
endif()

set( CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "--expt-relaxed-constexpr" )

set(RSA_MPI_INCLUDE_DIRS ${ONIKA_INCLUDE_DIRS} ${rsa_mpi_DIR}/include ${rsa_mpi_DIR}/prototype/include)

add_subdirectory(src)

add_executable(rsa_mpi rsa_mpi.cpp)
target_compile_options(rsa_mpi PUBLIC ${ONIKA_COMPILE_OPTIONS})
target_compile_definitions(rsa_mpi PUBLIC ${ONIKA_COMPILE_DEFINITIONS})
target_include_directories(rsa_mpi PUBLIC ${RSA_MPI_INCLUDE_DIRS})
target_link_directories(rsa_mpi PUBLIC ${ONIKA_LIBRARY_DIRS}) #${exaDEM_LIBRARY_DIRS})
target_link_libraries(rsa_mpi PUBLIC ${ONIKA_LIBRARIES} exanbCore)

set(ONIKA_PLUGIN_PATH "${CMAKE_CURRENT_BINARY_DIR}/plugins")
string(REPLACE ";" ":" ONIKA_PLUGIN_PATH "${ONIKA_PLUGIN_PATH}")

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/setup-env.sh
  "#!/bin/bash\n"
  "export onika_DIR=${onika_DIR}\n"
  "export exaNBody_DIR=${exaNBody_DIR}\n"
  "source ${onika_DIR}/bin/setup-env.sh\n"
  "export ONIKA_CONFIG_PATH=${CMAKE_CURRENT_SOURCE_DIR}/data/config\n"
  "export ONIKA_PLUGIN_PATH=${ONIKA_PLUGIN_PATH} \n")

add_subdirectory(test)

