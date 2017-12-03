# Define a variable basic_boundary_conditions_source_files containing
# a list of the source files for the basic_boundary_conditions library
set(basic_boundary_conditions_source_files
    BoundaryUtilityStrategy.cpp
    BasicCartesianBoundaryUtilities1.cpp
    BasicCartesianBoundaryUtilities2.cpp
    BasicCartesianBoundaryUtilities3.cpp
)

# Create a library called basic_boundary_conditions which includes the 
# source files defined in basic_boundary_conditions_source_files
add_library(basic_boundary_conditions ${basic_boundary_conditions_source_files})

TARGET_LINK_LIBRARIES(basic_boundary_conditions
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(basic_boundary_conditions z)
endif()
