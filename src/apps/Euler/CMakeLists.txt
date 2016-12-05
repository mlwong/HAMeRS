# Define a variable Euler_source_files containing a list of the source
# files for the Euler library
set(Euler_source_files
    Euler.cpp
    EulerBoundaryConditions.cpp
    EulerInitialConditions.cpp
)

# Create a library called Euler which includes the source files defined
# in Euler_source_files
add_library(Euler ${Euler_source_files})

TARGET_LINK_LIBRARIES(Euler
    convective_flux_reconstructors refinement_taggers flow_models
    Runge_Kutta_patch_strategy Runge_Kutta_level_integrator
    SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
    SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
    hdf5 gfortran)
