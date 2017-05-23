# Define a variable variable_fill_patterns_source_files containing
# a list of the source files for the variable_fill_patterns library
set(variable_fill_patterns_source_files
    LayerCellVariableFillPattern.cpp
    LayerCellNoCornersVariableFillPattern.cpp
    CellCornersVariableFillPattern.cpp
)

# Create a library called variable_fill_patterns which includes the 
# source files defined in variable_fill_patterns_source_files
add_library(variable_fill_patterns ${variable_fill_patterns_source_files})

TARGET_LINK_LIBRARIES(variable_fill_patterns
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
