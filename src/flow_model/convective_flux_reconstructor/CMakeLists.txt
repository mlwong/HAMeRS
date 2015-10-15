# Define a variable convective_flux_reconstructor_source_files containing
# a list of the source files for the convective_flux_reconstructor library
set(convective_flux_reconstructor_source_files
    ConvectiveFluxReconstructorFirstOrderHLLC.cpp
    ConvectiveFluxReconstructorLLF.cpp
    ConvectiveFluxReconstructorWENO-JS5-LLF.cpp
    ConvectiveFluxReconstructorWENO-CU6-M2-LLF.cpp
    ConvectiveFluxReconstructorWCNS-JS5-HLLC-HLL.cpp
    ConvectiveFluxReconstructorWCNS-CU6-M2-HLLC-HLL.cpp
    ConvectiveFluxReconstructorWCNS-HW6-HLLC-HLL.cpp
    ConvectiveFluxReconstructorWCNS-HW6-LD-HLLC-HLL.cpp
)

# Create a library called convective_flux_reconstructor which includes the 
# source files defined in convective_flux_reconstructor_source_files
add_library(convective_flux_reconstructor ${convective_flux_reconstructor_source_files})

TARGET_LINK_LIBRARIES(convective_flux_reconstructor
  Riemann_solver algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
