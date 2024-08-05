#!/bin/bash
#
# Environment definition for installation
#
# marc.josien@cea.fr
#
###############################################################################################


if [[ "${HOSTNAME}" == *"p00leiades"* ]]
then
    export PYBIND_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/copy_pybind 
else
    export PYBIND_REPO=https://github.com/pybind/pybind11
fi




