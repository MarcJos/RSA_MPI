set(PROTOTYPE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/../..)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../run_msp_test.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/run_msp_test.sh
  @ONLY)

add_test(
  NAME prototype_exanb_naif_3D
  COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/run_msp_test.sh ${CMAKE_CURRENT_SOURCE_DIR}/../../naif_3D_exanb_mpi.msp)

set_tests_properties(prototype_exanb_naif_3D PROPERTIES
  TIMEOUT 120
  PASS_REGULAR_EXPRESSION "iteration 10")
