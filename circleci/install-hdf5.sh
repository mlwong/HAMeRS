#!/bin/sh
set -e

wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz
tar xzf hdf5-1.14.6.tar.gz
cd hdf5-1.14.6
./configure --prefix=$HDF5_ROOT --without-zlib
make -j 4 > /dev/null 2>&1
make install
