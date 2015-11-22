# Define a variable Euler_boundary_conditions_source_files containing
# a list of the source files for the Euler_boundary_conditions library
set(Euler_boundary_conditions_source_files
    EulerBoundaryConditions.cpp
)

# Create a library called Euler_boundary_conditions which includes the 
# source files defined in Euler_boundary_conditions_source_files
add_library(Euler_boundary_conditions ${Euler_boundary_conditions_source_files})

TARGET_LINK_LIBRARIES(Euler_boundary_conditions
  equation_of_state algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
