cmake_minimum_required(VERSION 3.0)
project(rsa_mpi)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE "Release")

#packages
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# code coverage
option(CODE_COVERAGE OFF)
if(CODE_COVERAGE)
	set(GCC_COVERAGE_COMPILE_FLAGS "-g -O0 -coverage -fprofile-arcs -ftest-coverage")
	set(GCC_COVERAGE_LINK_FLAGS    "-coverage -lgcov")
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
	set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
endif()

add_subdirectory(src)

add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_executable(rsa command_line.cpp)
target_include_directories(rsa PUBLIC include/)

target_link_libraries(rsa PUBLIC MPI::MPI_CXX interfMatools OpenMP::OpenMP_CXX)

target_compile_options(rsa PRIVATE -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual -Wreturn-type -Wfloat-equal -Wno-endif-labels -Wmissing-format-attribute -Wno-multichar -Wno-deprecated-declarations -Wpacked -Wredundant-decls -Wunknown-pragmas -Wundef -Wreorder -Wno-conversion -Wdouble-promotion -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Wstrict-null-sentinel -Wtrampolines -Wvector-operation-performance -Wignored-qualifiers -Wsequence-point -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=format -Wswitch-bool -Wlogical-not-parentheses -Wsizeof-array-argument -Wbool-compare -Winline -Wshift-negative-value -Wshift-overflow -Wtautological-compare -Wnull-dereference -Wduplicated-cond -Wmisleading-indentation -Wduplicated-branches -Wregister -Wint-in-bool-context -Wsuggest-override -Wextra -pedantic -Wswitch -Wmissing-declarations -Wzero-as-null-pointer-constant -Wlong-long -Wno-switch-enum -Wno-covered-switch-default -Wno-covered-switch-default -Wno-sign-compare -Wno-missing-declarations -Wno-weak-vtables -Wno-unused-macros -Wsuggest-final-types -Wsuggest-final-methods)

add_subdirectory(Interface_python)

enable_testing()
add_subdirectory(test)

add_subdirectory(scripts)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION include
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/share/cmake/
	DESTINATION .
)
