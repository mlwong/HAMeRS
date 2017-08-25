#!/bin/sh

git clone --branch hdf5_1_8_14 https://git.hdfgroup.org/scm/hdffv/hdf5.git hdf5_1_8_14
cd hdf5_1_8_14
./configure --prefix=$HDF5_ROOT
make > /dev/null 2>&1
make install
