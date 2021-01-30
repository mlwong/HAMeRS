# HAMeRS: Hydrodynamics Adaptive Mesh Refinement Simulator #

[![Build Status](https://travis-ci.org/mlwong/HAMeRS.svg?branch=master)](https://travis-ci.org/mlwong/HAMeRS)

[HAMeRS](https://fpal.stanford.edu/hamers) is a compressible Navier-Stokes/Euler solver with the patch-based adaptive mesh refinement (AMR) technique. The parallelization of the code and all the construction, management and storage of cells are facilitated by the [Structured Adaptive Mesh Refinement Application Infrastructure](https://computation.llnl.gov/project/SAMRAI/) (SAMRAI) from the [Lawrence Livermore National Laboratory](https://www.llnl.gov/) (LLNL).

The code consists of various explicit high-order finite difference shock-capturing WCNS's (Weighted Compact Nonlinear Schemes) for capturing shock waves, material interfaces, and turbulent features. The AMR algorithm implemented is based on the one developed by Berger et al.

### How do I get set up? ###

Git is used for the version control of the code. To install Git on Debian-based distribution like Ubuntu, try apt-get:

```
sudo apt-get install git-all
```

To compile the code, in general all you need is to use [CMake](https://cmake.org/). For example, after cloning the repository with `git clone`:

```
cd HAMeRS
mkdir build
cd build
cmake ..
make
```

The compilers to be used to compile C, C++ and Fortran parts of HAMeRS can be chosen by setting the environment variables `CC`, `CXX` and `F77` respectively before running CMake. For example, to use the default MPI compilers, you can run:
```
export CC=mpicc
export CXX=mpicxx
export F77=mpif77
```

To run the code, you need to provide the input file:

```
src/exec/main <input filename>
```

To restart a simulation, you need to provide restart directory and restore number in addition to the input file:

```
src/exec/main <input filename> <restart dir> <restore number>
```

To run the code in parallel, you need MPI. You can try mpirun:

```
mpirun -np <number of processors> src/exec/main <input filename>
```

### What libraries do I need? ###

HAMeRS relies on [HDF5](https://support.hdfgroup.org/HDF5/), [Boost](https://www.boost.org/) and [SAMRAI](https://computation.llnl.gov/projects/samrai). Before installing HAMeRS, it is required to set up the environment variables for [CMake](https://cmake.org/) to look for the locations of the libraries.

To set up HDF5:
```
export HDF5_ROOT=<path to the directory of HDF5>
```

To set up Boost:
```
export BOOST_ROOT=<path to the directory of Boost>
```

To set up SAMRAI:
```
export SAMRAI_ROOT=<path to the directory of SAMRAI>
```

HAMeRS has already been successfully tested with HDF5-1.8, Boost-1.60 and SAMRAI-3.11.2.

Note that SAMRAI does not depend on the Boost library anymore since version 3.12.0. Please install HAMeRS without Boost library dependency using the CMake flag `-DHAMERS_USE_BOOST=OFF` when SAMRAI verison is equal to or greater than 3.12.0.

### How do I change the problem? ###

To change the problem that you want to run for an application, e.g. the Euler application, just simply link the corresponding initial conditions cpp symlink (`EulerInitialConditions.cpp` in `src/apps/Euler`) to the actual problem file using `ln -sf <absolute path to .cpp file containing problem's initial conditions> EulerInitialConditions.cpp`. If the problem has special boundary conditions, the user can supply the boundary conditions with `ln -sf <absolute path to .cpp file containing problem's user-coded boundary conditions> EulerSpecialBoundaryConditions.cpp`. There are some initial conditions and boundary conditions files from different example problems in the `problems` folder.

### Are there more tips and tutorials on how to compile and run the code?

Please have a look at the [Wiki page](https://github.com/mlwong/HAMeRS/wiki).

### Who do I talk to? ###

The code is managed by the previous PhD graduate Man-Long Wong (mlwong@alumni.stanford.edu) of the [Flow Physics and Aeroacoustics Laboratory](https://fpal.stanford.edu/) (FPAL)  at the [Department of Aeronautics and Astronautics](https://aa.stanford.edu/) of [Stanford University](https://www.stanford.edu/).

### Copyright ###
HAMeRS is licensed under a GNU Lesser General Public License v3.0.

If you find this work useful, please consider citing the author's dissertation:

    @phdthesis{wong2019thesis,
    title={High-order shock-capturing methods for study of shock-induced turbulent mixing with adaptive mesh refinement simulations},
    author={Wong, Man Long},
    year={2019},
    school={Stanford University}
    }
