# Define a variable equations_of_shear_viscosity_source_files containing
# a list of the source files for the equations_of_shear_viscosity library
set(equations_of_shear_viscosity_source_files
    EquationOfShearViscosityMixingRulesManager.cpp
    ./constant/EquationOfShearViscosityConstant.cpp
    ./constant/EquationOfShearViscosityMixingRulesConstant.cpp
    ./Chapman-Enskog/EquationOfShearViscosityChapmanEnskog.cpp
    ./Chapman-Enskog/EquationOfShearViscosityMixingRulesChapmanEnskog.cpp
)

# Create a library called equations_of_shear_viscosity which includes the 
# source files defined in equations_of_shear_viscosity_source_files
add_library(equations_of_shear_viscosity ${equations_of_shear_viscosity_source_files})

TARGET_LINK_LIBRARIES(equations_of_shear_viscosity algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
