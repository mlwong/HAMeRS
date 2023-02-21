#!/bin/sh

git clone --branch hdf5_1_8 https://github.com/HDFGroup/hdf5.git hdf5_1_8
cd hdf5_1_8
./configure --prefix=$HDF5_ROOT --without-zlib
make -j 4 > /dev/null 2>&1
make install
