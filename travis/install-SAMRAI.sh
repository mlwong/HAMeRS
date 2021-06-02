#!/bin/sh

cd ${HOME}
wget https://github.com/LLNL/SAMRAI/archive/refs/tags/v-4-0-2.tar.gz
gunzip v-4-0-2.tar.gz
mkdir SAMRAI-v4.0.2
cd SAMRAI-v4.0.2
tar xvf ../v-4-0-2.tar > /dev/null 2>&1
mkdir objs
cd objs
cmake ../SAMRAI-v-4-0-2 -DCMAKE_INSTALL_PREFIX=$SAMRAI_ROOT_NO_BOOST -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif77 -DENABLE_HDF5=ON -DHDF5_DIR=$HDF5_ROOT
make
make install
