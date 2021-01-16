# Define a variable flow_models_source_files containing
# a list of the source files for the flow_models library
set(flow_models_source_files
    FlowModelManager.cpp
    FlowModel.cpp
    FlowModelDiffusiveFluxUtilities.cpp
    FlowModelSourceUtilities.cpp
    ./MPI_helpers/FlowModelMPIHelperMaxMin.cpp
    ./MPI_helpers/FlowModelMPIHelperAverage.cpp
    ./MPI_helpers/FlowModelMPIHelperCorrelation.cpp
    ./single-species/FlowModelSingleSpecies.cpp
    ./single-species/FlowModelBasicUtilitiesSingleSpecies.cpp
    ./single-species/FlowModelBoundaryUtilitiesSingleSpecies.cpp
    ./single-species/FlowModelDiffusiveFluxUtilitiesSingleSpecies.cpp
    ./single-species/FlowModelRiemannSolverSingleSpecies.cpp
    ./single-species/FlowModelSourceUtilitiesSingleSpecies.cpp
    ./single-species/FlowModelStatisticsUtilitiesSingleSpecies.cpp
    ./single-species/Riemann_solvers/FlowModelRiemannSolverSingleSpeciesHLLC.cpp
    ./single-species/Riemann_solvers/FlowModelRiemannSolverSingleSpeciesHLLC-HLL.cpp
    ./four-eqn_conservative/FlowModelFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelBasicUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelBoundaryUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelDiffusiveFluxUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelRiemannSolverFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelSourceUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/Riemann_solvers/FlowModelRiemannSolverFourEqnConservativeHLLC.cpp
    ./four-eqn_conservative/Riemann_solvers/FlowModelRiemannSolverFourEqnConservativeHLLC-HLL.cpp
    ./five-eqn_Allaire/FlowModelFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelBasicUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelBoundaryUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelRiemannSolverFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelSourceUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelStatisticsUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/Riemann_solvers/FlowModelRiemannSolverFiveEqnAllaireHLLC.cpp
    ./five-eqn_Allaire/Riemann_solvers/FlowModelRiemannSolverFiveEqnAllaireHLLC-HLL.cpp
)

# Create a library called flow_models which includes the 
# source files defined in flow_models_source_files
add_library(flow_models ${flow_models_source_files})

TARGET_LINK_LIBRARIES(flow_models MPI_helpers
  visit_data_writer mixing_rules patch_hierarchies
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  ${HDF5_C_LIBRARIES} gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(flow_models z)
endif()
