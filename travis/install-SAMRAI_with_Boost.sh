#!/bin/sh

cd ${HOME}
wget https://github.com/LLNL/SAMRAI/archive/refs/tags/v-3-11-2.tar.gz
gunzip v-3-11-2.tar.gz
mkdir SAMRAI-v3.11.2
cd SAMRAI-v3.11.2
tar xvf ../v-3-11-2.tar > /dev/null 2>&1
mkdir objs
cd objs
sh ../SAMRAI-v-3-11-2/configure --prefix=$SAMRAI_ROOT_WITH_BOOST --enable-opt --with-CXX=mpicxx --with-CC=mpicc --with-F77=mpif77 --with-boost --with-hdf5=$HDF5_ROOT
make library
make install
