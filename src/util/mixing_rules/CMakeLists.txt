# Define a variable mixing_rules_source_files containing
# a list of the source files for the mixing_rules library
set(mixing_rules_source_files
    ./equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.cpp
    ./equations_of_shear_viscosity/constant/EquationOfShearViscosityConstant.cpp
    ./equations_of_shear_viscosity/constant/EquationOfShearViscosityMixingRulesConstant.cpp
    ./equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityChapmanEnskog.cpp
    ./equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityMixingRulesChapmanEnskog.cpp
    ./equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.cpp
    ./equations_of_bulk_viscosity/constant/EquationOfBulkViscosityConstant.cpp
    ./equations_of_bulk_viscosity/constant/EquationOfBulkViscosityMixingRulesConstant.cpp
    ./equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityCramer.cpp
    ./equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityMixingRulesCramer.cpp
    ./equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.cpp
    ./equations_of_thermal_conductivity/constant/EquationOfThermalConductivityConstant.cpp
    ./equations_of_thermal_conductivity/constant/EquationOfThermalConductivityMixingRulesConstant.cpp
    ./equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityPrandtl.cpp
    ./equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityMixingRulesPrandtl.cpp
    ./equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRulesManager.cpp
    ./equations_of_mass_diffusivity/constant/EquationOfMassDiffusivityMixingRulesConstant.cpp
    ./equations_of_mass_diffusivity/Reid/EquationOfMassDiffusivityMixingRulesReid.cpp
    )

# Create a library called mixing_rules which includes the 
# source files defined in mixing_rules_source_files
add_library(mixing_rules ${mixing_rules_source_files})

TARGET_LINK_LIBRARIES(mixing_rules equations_of_state
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(mixing_rules z)
endif()
