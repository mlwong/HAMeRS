#!/bin/sh

git clone --branch hdf5_1_8 https://git.hdfgroup.org/scm/hdffv/hdf5.git hdf5_1_8
cd hdf5_1_8
./configure --prefix=$HDF5_ROOT --without-zlib
make > /dev/null 2>&1
make install
