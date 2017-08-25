# HAMeRS: Hydrodynamics Adaptive Mesh Refinement Simulator #

[![Build Status](https://travis-ci.com/mlwong/HAMeRS.svg?token=qhzno8A188sa9Lwm5JDg&branch=master)](https://travis-ci.com/mlwong/HAMeRS)

[HAMeRS](https://fpal.stanford.edu/hamers) is a compressible Navier-Stokes/Euler solver with the block-structured adaptive mesh refinement (AMR) technique. The parallelization of the code and all the construction, management and storage of cells are facilitated by the [Structured Adaptive Mesh Refinement Application Infrastructure](https://computation.llnl.gov/project/SAMRAI/) (SAMRAI) from the [Lawrence Livermore National Laboratory](https://www.llnl.gov/) (LLNL).

The code consists of the families of explicit high-order finite difference shock-capturing WCNS's (Weighted Compact Nonlinear Schemes) for capturing shock waves, material interfaces, and turbulent features. The AMR algorithm implemented is based on the one developed by Berger et al.

### How do I get set up? ###

Git is used for the version control of the code. To install Git on Debian-based distribution like Ubuntu, try apt-get:

```
sudo apt-get install git-all
```

To compile the code, in general all you need is to use cmake and then make. For example:

```
mkdir build
cd build
cmake ..
make
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

### Who do I talk to? ###

The code is managed by Man-Long Wong (wongml@stanford.edu) of the [Flow Physics and Aeroacoustics Laboratory](https://fpal.stanford.edu/) (FPAL)  at the [Department of Aeronautics and Astronautics](https://aa.stanford.edu/) of [Stanford University](https://www.stanford.edu/).

### Copyright ###
HAMeRS is licensed under a  GNU Lesser General Public License.
