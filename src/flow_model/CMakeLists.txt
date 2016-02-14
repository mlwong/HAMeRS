# Define a variable flow_model_manager_source_files containing
# a list of the source files for the flow_model_manager library
set(flow_model_manager_source_files
    FlowModelManager.cpp
)

# Create a library called flow_model_manager which includes the 
# source files defined in flow_model_manager_source_files
add_library(flow_model_manager ${flow_model_manager_source_files})

TARGET_LINK_LIBRARIES(flow_model_manager initial_conditions
  Euler_boundary_conditions convective_flux_reconstructor
  refinement_taggers equation_of_state algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
