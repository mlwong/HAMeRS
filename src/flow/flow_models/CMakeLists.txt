# Define a variable flow_models_source_files containing
# a list of the source files for the flow_models library
set(flow_models_source_files
    FlowModelManager.cpp
    FlowModelSingleSpecies.cpp
    FlowModelFourEqnConservative.cpp
    FlowModelFiveEqnAllaire.cpp
)

# Create a library called flow_models which includes the 
# source files defined in flow_models_source_files
add_library(flow_models ${flow_models_source_files})

TARGET_LINK_LIBRARIES(flow_models Riemann_solvers
  equations_of_state algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
