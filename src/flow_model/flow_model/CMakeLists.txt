# Define a variable flow_model_source_files containing
# a list of the source files for the flow_model library
set(flow_model_source_files
    FlowModelSingleSpecies.cpp
    FlowModelFourEqnConservative.cpp
    FlowModelFiveEqnAllaire.cpp
)

# Create a library called flow_model which includes the 
# source files defined in flow_model_source_files
add_library(flow_model ${flow_model_source_files})

TARGET_LINK_LIBRARIES(flow_model Riemann_solver algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
