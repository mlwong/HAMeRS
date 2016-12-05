# Define a variable diffusive_flux_reconstructors_source_files containing
# a list of the source files for the diffusive_flux_reconstructors library
set(diffusive_flux_reconstructors_source_files
    DiffusiveFluxReconstructorManager.cpp
    ./sixth_order/DiffusiveFluxReconstructorSixthOrder.cpp
)

# Create a library called diffusive_flux_reconstructors which includes the 
# source files defined in diffusive_flux_reconstructors_source_files
add_library(diffusive_flux_reconstructors ${diffusive_flux_reconstructors_source_files})

TARGET_LINK_LIBRARIES(diffusive_flux_reconstructors
  flow_models algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
