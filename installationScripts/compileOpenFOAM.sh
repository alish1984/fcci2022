#!/usr/bin/env bash

TARGET_DIR="${HOME}/OpenFOAM"

[[ ! -d ${TARGET_DIR} ]] && mkdir -p $TARGET_DIR

cd ${TARGET_DIR} \
&& wget https://dl.openfoam.com/source/v2106/OpenFOAM-v2106.tgz \
&& wget https://dl.openfoam.com/source/v2106/ThirdParty-v2106.tgz \
&& tar -xzf OpenFOAM-v2106.tgz \
&& tar -xzf ThirdParty-v2106.tgz \
&& source ${TARGET_DIR}/OpenFOAM-v2106/etc/bashrc \
&& cp $WM_PROJECT_DIR/etc/config.sh/CGAL $WM_PROJECT_DIR/etc/config.sh/CGAL.back \
&& sed -i s/boost_version=boost_1_66_0/boost_version=boost_1_70_0/g $WM_PROJECT_DIR/etc/config.sh/CGAL \
&& sed -i '49 i GCC_VERS1=${GCC_VERS1:-"${COMPILER}`gcc --version  | awk '\''/gcc/{print $4}'\''`"}' $WM_PROJECT_DIR/etc/config.sh/CGAL \
&& sed -e '/export BOOST_ARCH_PATH=/s/^/#/g' -i $WM_PROJECT_DIR/etc/config.sh/CGAL \
&& sed -i '50 i export BOOST_ARCH_PATH=${HOME}/NumericalLibraries/boost/boost_1_70_0-gcc-${GCC_VERS1}' $WM_PROJECT_DIR/etc/config.sh/CGAL \
&& foamSystemCheck >log.foamSystemCheck 2>&1 \
&& cd ${TARGET_DIR}/OpenFOAM-v2106 \
&& ./Allwmake -j -s -q -l -k \
&& ./Allwmake -s -q > log.build2 2>&1 \
&& foamInstallationTest > log.foamInstallationTest 2>&1 \
&& mkdir -p $FOAM_RUN \
&& cd $FOAM_RUN \
&& cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily ./ \
&& cd pitzDaily \
&& blockMesh \
&& simpleFoam

echo "alias ofv2106='source ${TARGET_DIR}/OpenFOAM-v2106/etc/bashrc'" >> ${HOME}/.bashrc

echo "---------------- - ------------------ - ----------------"
echo "----------------  compilation finished  ----------------"
echo "-------------  Please check the log files  -------------"
echo "---- to make sure all targets updated successfully  ----"
echo "---------------- - ------------------ - ----------------"


