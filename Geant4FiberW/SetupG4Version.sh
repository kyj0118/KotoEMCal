#!/bin/bash

# Geant4 version setting
#g4version=10.6.3
g4version=10.7.1
source /sw/packages/geant4/$g4version/bin/geant4.sh
source /sw/packages/geant4/$g4version/share/Geant4-${g4version}/geant4make/geant4make.sh

# cmake3
export CMAKE=/gpfs/fs03/had/koto/ps3/klea/work/yjkim/public/tools/cmake3
export PATH=$CMAKE/bin:$PATH
export LD_LIBRARY_PATH=$GCC/lib64:$LD_LIBRARY_PATH

# Compiler version setting
module load gcc/930
export GCC=/opt/gcc-9.3.0
export CC=$GCC/bin/gcc
export LD_LIBRARY_PATH=$GCC/lib64:$LD_LIBRARY_PATH
export PATH=$GCC/bin:$PATH

# Root version setting
export ROOTSYS=/sw/packages/root/pro
source $ROOTSYS/bin/thisroot.sh
