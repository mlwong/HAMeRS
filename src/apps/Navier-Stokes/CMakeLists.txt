# Define a variable Navier_Stokes_source_files containing a list of the source
# files for the Navier_Stokes library
set(Navier_Stokes_source_files
    NavierStokes.cpp
    NavierStokesBoundaryConditions.cpp
    NavierStokesInitialConditions.cpp
)

# Create a library called Navier_Stokes which includes the source files defined
# in Navier_Stokes_source_files
add_library(Navier_Stokes ${Navier_Stokes_source_files})

TARGET_LINK_LIBRARIES(Navier_Stokes
    convective_flux_reconstructors diffusive_flux_reconstructors
    nonconservative_diffusive_flux_divergence_operators refinement_taggers
    flow_models Runge_Kutta_patch_strategy Runge_Kutta_level_integrator
    SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
    SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
    ${HDF5_C_LIBRARIES} gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(Navier_Stokes z)
endif()
