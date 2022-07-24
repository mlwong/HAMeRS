#!/bin/sh

cd ${HOME}
wget https://github.com/LLNL/SAMRAI/releases/download/v-4-1-0/SAMRAI-v4.1.0.tar.gz
gunzip SAMRAI-v4.1.0.tar.gz
mkdir SAMRAI-v4.1.0
cd SAMRAI-v4.1.0
tar xvf ../SAMRAI-v4.1.0.tar > /dev/null 2>&1
mkdir objs
cd objs
cmake ../SAMRAI -DCMAKE_INSTALL_PREFIX=$SAMRAI_ROOT_NO_BOOST -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif77 -DENABLE_HDF5=ON -DHDF5_DIR=$HDF5_ROOT
make -j 4
make install
