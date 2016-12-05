# Define a variable utils_source_files containing
# a list of the source files for the utils library
set(utils_source_files
    ./gradient_sensors/GradientSensorJameson.cpp
    ./wavelet_transform/WaveletTransformHarten.cpp
)

# Create a library called utils which includes the 
# source files defined in utils_source_files
add_library(utils ${utils_source_files})

TARGET_LINK_LIBRARIES(utils mixing_rules
  basic_boundary_conditions algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
