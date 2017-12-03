# Define a variable Runge_Kutta_patch_strategy_source_files containing 
# a list of the source files for the Runge_Kutta_patch_strategy library
set(Runge_Kutta_patch_strategy_source_files
    RungeKuttaPatchStrategy.cpp
)

# Create a library called Runge_Kutta_patch_strategy which includes the 
# source files defined in Runge_Kutta_patch_strategy_source_files
add_library(Runge_Kutta_patch_strategy ${Runge_Kutta_patch_strategy_source_files})

TARGET_LINK_LIBRARIES(Runge_Kutta_patch_strategy 
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  ${HDF5_C_LIBRARIES} gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(Runge_Kutta_patch_strategy z)
endif()
