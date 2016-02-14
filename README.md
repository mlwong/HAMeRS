# README #

### What is this repository for? ###

HOAMR (High-Order Adaptive Mesh Refinement) is a compressible Navier-Stokes/Euler solver with the AMR (Adaptive Mesh Refinement) technique. The parallelization of the code and all the construction, management and storage of cells are facilitated by the SAMRAI (Structured Adaptive Mesh Refinement Application Infrastructure) library from the Lawrence Livermore National Lab (computation.llnl.gov/project/SAMRAI/).

The code consists of a choice of explicit WCNS (Weighted Compact Nonlinear Scheme) and WENO (Weighted Essentially Non-Oscillatory) high-order shock-capturing schemes for simulating shock waves, material interfaces, and turbulent features. The AMR algorithm implemented is based on the algorithm developed by Berger et al. The figure below summarizes the algorithm:

![AMR_flowchart.png](https://bitbucket.org/repo/zzaMX8/images/2144227433-AMR_flowchart.png)


### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions


### Who do I talk to? ###

The code is managed by Man Long, Wong of the UFPA lab (Unsteady Flow Physics and Aeroacoustics Laboratory) at Stanford University. If you have any questions, please feel free to contact Man Long at wongml@stanford.edu