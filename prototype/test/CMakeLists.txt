# licensed to the apache software foundation (asf) under one
# or more contributor license agreements.  see the notice file
# distributed with this work for additional information
# regarding copyright ownership.  the asf licenses this file
# to you under the apache license, version 2.0 (the
# "license"); you may not use this file except in compliance
# with the license.  you may obtain a copy of the License at
#
#   http://www.apache.org/licenses/license-2.0
#
# unless required by applicable law or agreed to in writing,
# software distributed under the license is distributed on an
# "as is" basis, without warranties or conditions of any
# kind, either express or implied.  see the license for the
# specific language governing permissions and limitations
# under the license.

# runs a built rsa_mpi + its plugins against a .msp scenario, through the
# same setup-env.sh (ONIKA_CONFIG_PATH / ONIKA_PLUGIN_PATH / LD_LIBRARY_PATH)
# used for manual runs, so ctest exercises the exact same onika app.
set(PROTOTYPE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/..)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/run_msp_test.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/run_msp_test.sh
  @ONLY)

add_test(
  NAME prototype_naif_3D
  COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/run_msp_test.sh ${CMAKE_CURRENT_SOURCE_DIR}/../naif_3D.msp)

# the naive RSA loop is meant to terminate on its own (either n_draw is
# reached or the radius generator runs out of radii to place); a hang here
# is a regression in the loop-continuation logic, so fail fast instead of
# stalling the whole test suite.
set_tests_properties(prototype_naif_3D PROPERTIES
  TIMEOUT 120
  PASS_REGULAR_EXPRESSION "Number of Spheres")
