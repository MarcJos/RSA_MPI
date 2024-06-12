#!/bin/bash
#
# Script for downloading the prerequisite
#
# Should be executed in the root of the directory
#
# marc.josien@cea.fr
#
###############################################################################################

source Installation/Install_environment.sh

## get  pybind11
cd Interface_python
rm -rf pybind11
git clone $PYBIND_REPO pybind11
cd ../

