configure_file("test_py.sh.in" ${CMAKE_BINARY_DIR}/test/test_py.sh @ONLY)    
add_test(NAME test_python COMMAND bash test_py.sh)

# Notre exécutable
add_executable(test_helper_setter test_helper_setter.cpp)
target_include_directories(test_helper_setter PUBLIC ../include/)
add_test(NAME test_helper_setter_dim_1 COMMAND ./test_helper_setter 1)
add_test(NAME test_helper_setter_dim_2 COMMAND ./test_helper_setter 2)
add_test(NAME test_helper_setter_dim_3 COMMAND ./test_helper_setter 3)
add_test(NAME test_helper_setter_dim_10 COMMAND ./test_helper_setter 10)

add_executable(test_data_storage test_data_storage.cpp)
target_link_libraries(test_data_storage MPI::MPI_CXX )
target_include_directories(test_data_storage PUBLIC ../include/)
add_test(NAME fill_data_storage_serial COMMAND ./test_data_storage)

add_executable(test_minimal_decoration test_decoration.cpp)
target_link_libraries(test_minimal_decoration MPI::MPI_CXX )
target_include_directories(test_minimal_decoration PRIVATE ../include/)
add_test(NAME test_minimal_decoration_mpi COMMAND mpirun -n 4 ./test_minimal_decoration)

add_executable(test_3D_full test_3D_full.cpp)
target_link_libraries(test_3D_full MPI::MPI_CXX )
target_include_directories(test_3D_full PUBLIC ../include/)
add_test(NAME test_3D_full_serial COMMAND ./test_3D_full)
add_test(NAME test_3D_full_mpi COMMAND mpirun -np 4 ./test_3D_full)

add_executable(test_2D_full test_2D_full.cpp)
target_link_libraries(test_2D_full MPI::MPI_CXX )
target_include_directories(test_2D_full PUBLIC ../include/)
add_test(NAME test_2D_full_serial COMMAND ./test_2D_full)
add_test(NAME test_2D_full_mpi COMMAND mpirun -np 4 ./test_2D_full)

add_executable(test_3D_full_voxel test_3D_full_voxel.cpp)
target_link_libraries(test_3D_full_voxel MPI::MPI_CXX )
target_include_directories(test_3D_full_voxel PUBLIC ../include/)
add_test(NAME test_3D_full_voxel_serial COMMAND ./test_3D_full_voxel)
add_test(NAME test_3D_full_voxel_mpi COMMAND mpirun -np 4 ./test_3D_full_voxel)

add_executable(test_3D_radius test_3D_radius.cpp)
target_link_libraries(test_3D_radius MPI::MPI_CXX )
target_include_directories(test_3D_radius PUBLIC ../include/)
add_test(NAME test_3D_radius COMMAND mpirun -np 4 ./test_3D_radius)

add_executable(test_4D_full_voxel test_4D_full_voxel.cpp)
target_link_libraries(test_4D_full_voxel MPI::MPI_CXX )
target_include_directories(test_4D_full_voxel PUBLIC ../include/)
add_test(NAME test_4D_full_voxel_serial COMMAND ./test_4D_full_voxel)
add_test(NAME test_4D_full_voxel_mpi COMMAND mpirun -np 4 ./test_4D_full_voxel)

add_executable(test_2D_full_voxel test_2D_full_voxel.cpp)
target_link_libraries(test_2D_full_voxel MPI::MPI_CXX )
target_include_directories(test_2D_full_voxel PUBLIC ../include/)
add_test(NAME test_2D_full_voxel_serial COMMAND ./test_2D_full_voxel)
add_test(NAME test_2D_full_voxel_mpi COMMAND mpirun -np 4 ./test_2D_full_voxel)


add_executable(test_plugin_write_xyz test_plugin_write_xyz.cpp)
target_link_libraries(test_plugin_write_xyz MPI::MPI_CXX OpenMP::OpenMP_CXX )
target_include_directories(test_plugin_write_xyz PUBLIC ../include/ ../plugin)
add_test(NAME test_plugin_write_xyz COMMAND ./test_plugin_write_xyz 0.5 10 10 10)
add_test(NAME test_plugin_write_xyz_mpi COMMAND mpirun -np 4 ./test_plugin_write_xyz 0.5 10 10 10)

add_executable(test_split_concatenate test_split_concatenate.cpp)
target_include_directories(test_split_concatenate PUBLIC ../include/)
add_test(NAME test_split_concatenate COMMAND ./test_split_concatenate)

add_executable(test_cpp test_python.cpp)
target_link_libraries(test_cpp MPI::MPI_CXX)
target_include_directories(test_cpp PUBLIC ../include/)
add_test(NAME test_cpp COMMAND ./test_cpp)

add_test(NAME test_interface_2 COMMAND ./../rsa --dim 4 --radius 0.1 --seed 0 --size 10000 --inf 0, 0, 0, 0 --sup 1, 1, 1, 1)
