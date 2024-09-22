# Define a variable hyperviscosity_operators_source_files containing
# a list of the source files for the hyperviscosity_operators library
set(hyperviscosity_operators_source_files
    HyperviscosityOperator.cpp
)

# Create a library called hyperviscosity_operators which includes the 
# source files defined in hyperviscosity_operators_source_files
add_library(hyperviscosity_operators
  ${hyperviscosity_operators_source_files})

TARGET_LINK_LIBRARIES(hyperviscosity_operators
  flow_models
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  ${HDF5_C_LIBRARIES} gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(hyperviscosity_operators z)
endif()
