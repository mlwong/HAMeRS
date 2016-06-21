# Define a variable Riemann_solvers_source_files containing
# a list of the source files for the Riemann_solvers library
set(Riemann_solvers_source_files
    RiemannSolverSingleSpeciesHLLC.cpp
    RiemannSolverSingleSpeciesHLLC_HLL.cpp
    RiemannSolverFourEqnConservativeHLLC.cpp
    RiemannSolverFourEqnConservativeHLLC_HLL.cpp
    RiemannSolverFiveEqnAllaireHLLC.cpp
    RiemannSolverFiveEqnAllaireHLLC_HLL.cpp
)

# Create a library called Riemann_solvers which includes the 
# source files defined in Riemann_solvers_source_files
add_library(Riemann_solvers ${Riemann_solvers_source_files})

TARGET_LINK_LIBRARIES(Riemann_solvers
  equations_of_state algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
