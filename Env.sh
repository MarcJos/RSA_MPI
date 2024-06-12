#!/bin/bash
# 1) Get the path of this script (assumed to be in the root)
pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
if ([ -h "${SCRIPT_PATH}" ]); then
  while([ -h "${SCRIPT_PATH}" ]); do cd `dirname "$SCRIPT_PATH"`;
  SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null

# get the environment compiler of the installation
source ${SCRIPT_PATH}/Installation/Install_environment.sh
# Python library
export PYTHONPATH=${SCRIPT_PATH}/build/Interface_python:$PYTHONPATH
# build
export PATH=${PATH}:${SCRIPT_PATH}/build:${SCRIPT_PATH}/build/scripts
