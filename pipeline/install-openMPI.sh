#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>

set -e
set -x

os=`uname`

case "$os" in
    Darwin)
        brew update
        brew upgrade cmake
        brew upgrade gcc
		brew install openmpi
		;;

    Linux)
        sudo apt-get update -q
        sudo apt-get install -y gfortran
        wget --no-check-certificate https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.2.tar.gz
        tar -zxf openmpi-1.10.2.tar.gz
        cd openmpi-1.10.2
        sh ./configure --disable-dlopen --prefix=${HOME}/openmpi > /dev/null
        make install > /dev/null
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac