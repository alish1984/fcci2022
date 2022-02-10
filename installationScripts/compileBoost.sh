#!/usr/bin/env bash

TARGET_DIR="${HOME}/NumericalLibraries"
GCC_VERS=${GCC_VERS:-"${COMPILER}`gcc --version  | awk '/gcc/{print $4}'`"}

[[ ! -d ${TARGET_DIR} ]] && mkdir -p $TARGET_DIR

cd ${TARGET_DIR} \
&& mkdir -p ${TARGET_DIR}/boost/boost_1_70_0-gcc-${GCC_VERS} \
&& wget https://boostorg.jfrog.io/artifactory/main/release/1.70.0/source/boost_1_70_0.tar.gz\
&& tar -xzf boost_1_70_0.tar.gz -C ${TARGET_DIR}/boost/boost_1_70_0-gcc-${GCC_VERS} \
&& cd ${TARGET_DIR}/boost/boost_1_70_0-gcc-${GCC_VERS} \
&& mv boost_1_70_0 src && cd src && ls \
&& ./bootstrap.sh \
&& ./b2 install --prefix=${TARGET_DIR}/boost/boost_1_70_0-gcc-${GCC_VERS} 2>&1 | tee compile_boost1.70.0-gcc-${GCC_VERS}.log

echo "export LD_LIBRARY_PATH=${HOME}/NumericalLibraries/boost/boost_1_70_0-gcc-${GCC_VERS}/lib/:\$LD_LIBRARY_PATH" >> ${HOME}/.bashrc

echo "---------------- - ------------------ - ----------------"
echo "----------------  compilation finished  ----------------"
echo "-------------  Please check the log file  --------------"
echo "---- to make sure all targets updated successfully  ----"
echo "---------------- - ------------------ - ----------------"
