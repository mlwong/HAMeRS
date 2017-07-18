# Define a variable flow_models_source_files containing
# a list of the source files for the flow_models library
set(flow_models_source_files
    FlowModelManager.cpp
    FlowModel.cpp
    ./single-species/FlowModelSingleSpecies.cpp
    ./single-species/FlowModelBoundaryUtilitiesSingleSpecies.cpp
    ./single-species/FlowModelStatisticsUtilitiesSingleSpecies.cpp
    ./single-species/Riemann_solvers/RiemannSolverSingleSpeciesHLLC.cpp
    ./single-species/Riemann_solvers/RiemannSolverSingleSpeciesHLLC_HLL.cpp
    ./four-eqn_conservative/FlowModelFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelBoundaryUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.cpp
    ./four-eqn_conservative/Riemann_solvers/RiemannSolverFourEqnConservativeHLLC.cpp
    ./four-eqn_conservative/Riemann_solvers/RiemannSolverFourEqnConservativeHLLC_HLL.cpp
    ./five-eqn_Allaire/FlowModelFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelBoundaryUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/FlowModelStatisticsUtilitiesFiveEqnAllaire.cpp
    ./five-eqn_Allaire/Riemann_solvers/RiemannSolverFiveEqnAllaireHLLC.cpp
    ./five-eqn_Allaire/Riemann_solvers/RiemannSolverFiveEqnAllaireHLLC_HLL.cpp
)

# Create a library called flow_models which includes the 
# source files defined in flow_models_source_files
add_library(flow_models ${flow_models_source_files})

TARGET_LINK_LIBRARIES(flow_models mixing_rules patch_hierarchies 
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
