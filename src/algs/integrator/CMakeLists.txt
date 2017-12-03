# Define a variable Runge_Kutta_level_integrator_source_files containing
# a list of the source files for the Runge_Kutta_level_integrator library
set(Runge_Kutta_level_integrator_source_files
    ExtendedTagAndInitialize.cpp
    ExtendedTagAndInitializeConnectorWidthRequestor.cpp
    ExtendedTagAndInitStrategy.cpp
    RungeKuttaLevelIntegrator.cpp
    ./fortran/algs_upfluxsum3d.m4
    ./fortran/algs_upfluxsum3d.f
    ./fortran/algs_upfluxsum2d.m4
    ./fortran/algs_upfluxsum2d.f
    ./fortran/algs_upfluxsum1d.m4
    ./fortran/algs_upfluxsum1d.f
    ./fortran/algs_nodeouternodeops3d.m4
    ./fortran/algs_nodeouternodeops3d.f
    ./fortran/algs_nodeouternodeops2d.m4
    ./fortran/algs_nodeouternodeops2d.f
)

# Create a library called Runge_Kutta_level_integrator which includes the 
# source files defined in Runge_Kutta_level_integrator_source_files
add_library(Runge_Kutta_level_integrator ${Runge_Kutta_level_integrator_source_files})

TARGET_LINK_LIBRARIES(Runge_Kutta_level_integrator 
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(Runge_Kutta_level_integrator z)
endif()
