#!/bin/sh

cd ${HOME}
wget https://computation.llnl.gov/projects/samrai/download/SAMRAI-v3.11.2.tar.gz
gunzip SAMRAI-v3.11.2.tar.gz
mkdir SAMRAI-v3.11.2
cd SAMRAI-v3.11.2
tar xvf ../SAMRAI-v3.11.2.tar
mkdir objs
cd objs
sh ../SAMRAI/configure --prefix=$SAMRAI_ROOT --quiet --enable-opt --with-CXX=mpicxx --with-CC=mpicc --with-F77=mpif77 --with-boost --with-hdf5=${HOME}/hdf5
make library > /dev/null
make install > /dev/null
