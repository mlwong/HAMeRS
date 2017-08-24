#!/bin/sh

set -e
set -x

cd ${HOME}
wget https://computation.llnl.gov/projects/samrai/download/SAMRAI-v3.11.2.tar.gz
tar zxvf SAMRAI-v3.11.2.tar.gz
cd SAMRAI
mkdir objs
cd objs
sh ../SAMRAI/configure --enable-opt --with-CXX=mpicxx --with-CC=mpicc --with-F77=mpif77 --with-boost --with-hdf5