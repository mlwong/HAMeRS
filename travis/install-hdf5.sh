#!/bin/sh

git clone --branch hdf5_1_8_14 https://git.hdfgroup.org/scm/hdffv/hdf5.git hdf5_1_8_14
cd hdf5_1_8_14
./configure --prefix=${HOME}/hdf5 --quiet
make > /dev/null
make install > /dev/null
