# Define a variable Riemann_solver_source_files containing
# a list of the source files for the Riemann_solver library
set(Riemann_solver_source_files
    RiemannSolverHLLC.cpp
    RiemannSolverHLLC_HLL.cpp
    RiemannSolverSingleSpeciesHLLC.cpp
    RiemannSolverSingleSpeciesHLLC_HLL.cpp
    RiemannSolverFourEqnConservativeHLLC.cpp
    RiemannSolverFourEqnConservativeHLLC_HLL.cpp
    RiemannSolverFourEqnShyueHLLC.cpp
    RiemannSolverFourEqnShyueHLLC_HLL.cpp
    RiemannSolverFiveEqnAllaireHLLC.cpp
    RiemannSolverFiveEqnAllaireHLLC_HLL.cpp
)

# Create a library called Riemann_solver which includes the 
# source files defined in Riemann_solver_source_files
add_library(Riemann_solver ${Riemann_solver_source_files})

TARGET_LINK_LIBRARIES(Riemann_solver
  equation_of_state algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
