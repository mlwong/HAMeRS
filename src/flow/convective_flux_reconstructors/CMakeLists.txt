# Define a variable convective_flux_reconstructors_source_files containing
# a list of the source files for the convective_flux_reconstructors library
set(convective_flux_reconstructors_source_files
    ConvectiveFluxReconstructorManager.cpp
    ./first_order/ConvectiveFluxReconstructorFirstOrderLLF.cpp
    ./first_order/ConvectiveFluxReconstructorFirstOrderHLLC.cpp
    ./WCNS56/ConvectiveFluxReconstructorWCNS56-HLLC-HLL.cpp
    ./WCNS56/ConvectiveFluxReconstructorWCNS5-JS-HLLC-HLL.cpp
    ./WCNS56/ConvectiveFluxReconstructorWCNS5-Z-HLLC-HLL.cpp
    ./WCNS56/ConvectiveFluxReconstructorWCNS6-CU-M2-HLLC-HLL.cpp
    ./WCNS56/ConvectiveFluxReconstructorWCNS6-LD-HLLC-HLL.cpp
    ./WCNS56/ConvectiveFluxReconstructorWCNS6-Test.cpp
)

# Create a library called convective_flux_reconstructors which includes the 
# source files defined in convective_flux_reconstructors_source_files
add_library(convective_flux_reconstructors ${convective_flux_reconstructors_source_files})

TARGET_LINK_LIBRARIES(convective_flux_reconstructors flow_models
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(convective_flux_reconstructors z)
endif()
