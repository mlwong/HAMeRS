#!/usr/bin/env python3

import os, errno, pathlib
import argparse

# Euler application.
path_EulerInitialConditions_cpp         = ""
path_EulerSpecialBoundaryConditions_cpp = ""
path_EulerErrorStatistics_cpp           = ""

# Navier-Stokes application.
path_NavierStokesInitialConditions_cpp         = ""
path_NavierStokesSpecialBoundaryConditions_cpp = ""
path_NavierStokesErrorStatistics_cpp           = ""

# Flow models.
path_FlowModelSpecialSourceTerms_cpp                     = ""
path_FlowModelStatisticsUtilitiesSingleSpecies_cpp       = ""
path_FlowModelStatisticsUtilitiesFourEqnConservative_cpp = ""
path_FlowModelStatisticsUtilitiesFiveEqnAllaire_cpp      = ""

# Immersed boundaries.
path_ImmersedBoundaries_cpp = ""

### Default paths.

# Euler application.
path_EulerInitialConditions_cpp_default = "EulerInitialConditionsDefault.cpp"
path_EulerSpecialBoundaryConditions_cpp_default = "EulerSpecialBoundaryConditionsDefault.cpp"
path_EulerErrorStatistics_cpp_default = "EulerErrorStatisticsDefault.cpp"

# Navier-Stokes application.
path_NavierStokesInitialConditions_cpp_default = "NavierStokesInitialConditionsDefault.cpp"
path_NavierStokesSpecialBoundaryConditions_cpp_default = "NavierStokesSpecialBoundaryConditionsDefault.cpp"
path_NavierStokesErrorStatistics_cpp_default = "NavierStokesErrorStatisticsDefault.cpp"

# Flow models.
path_FlowModelSpecialSourceTerms_cpp_default = "FlowModelSpecialSourceTermsDefault.cpp"
path_FlowModelStatisticsUtilitiesSingleSpecies_cpp_default = "FlowModelStatisticsUtilitiesSingleSpeciesDefault.cpp"
path_FlowModelStatisticsUtilitiesFourEqnConservative_cpp_default = "FlowModelStatisticsUtilitiesFourEqnConservativeDefault.cpp"
path_FlowModelStatisticsUtilitiesFiveEqnAllaire_cpp_default = "FlowModelStatisticsUtilitiesFiveEqnAllaireDefault.cpp"

# Immersed boundaries.
path_ImmersedBoundaries_cpp_default = "ImmersedBoundariesDefault.cpp"

### Helper functions.
def symlinkForce(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)

def symlinkHAMeRSFile(path_src, path_src_default, path_des):
    if path_src == "":
        path = path_src_default
    else:
        path = path_src
    if path != os.readlink(path_des):
        symlinkForce(path, path_des)
        pathlib.Path(path_des)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reset", type=bool, help="enter true or false to reset all links to the default links",
                        nargs='?', default=False)
    args = parser.parse_args()
    if args.reset == True:
        # Euler application.
        path_EulerInitialConditions_cpp         = ""
        path_EulerSpecialBoundaryConditions_cpp = ""
        path_EulerErrorStatistics_cpp           = ""

        # Navier-Stokes application.
        path_NavierStokesInitialConditions_cpp         = ""
        path_NavierStokesSpecialBoundaryConditions_cpp = ""
        path_NavierStokesErrorStatistics_cpp           = ""

        # Flow models.
        path_FlowModelSpecialSourceTerms_cpp                     = ""
        path_FlowModelStatisticsUtilitiesSingleSpecies_cpp       = ""
        path_FlowModelStatisticsUtilitiesFourEqnConservative_cpp = ""
        path_FlowModelStatisticsUtilitiesFiveEqnAllaire_cpp      = ""

        # Immersed boundaries.
        path_ImmersedBoundaries_cpp = ""

    # Link files related to Euler applications.
    symlinkHAMeRSFile(path_EulerInitialConditions_cpp, path_EulerInitialConditions_cpp_default, "src/apps/Euler/EulerInitialConditions.cpp")
    symlinkHAMeRSFile(path_EulerSpecialBoundaryConditions_cpp, path_EulerSpecialBoundaryConditions_cpp_default, "src/apps/Euler/EulerSpecialBoundaryConditions.cpp")
    symlinkHAMeRSFile(path_EulerErrorStatistics_cpp, path_EulerErrorStatistics_cpp_default, "src/apps/Euler/EulerErrorStatistics.cpp")

    # Link files related to Navier-Stokes applications.
    symlinkHAMeRSFile(path_NavierStokesInitialConditions_cpp, path_NavierStokesInitialConditions_cpp_default, "src/apps/Navier-Stokes/NavierStokesInitialConditions.cpp")
    symlinkHAMeRSFile(path_NavierStokesSpecialBoundaryConditions_cpp, path_NavierStokesSpecialBoundaryConditions_cpp_default, "src/apps/Navier-Stokes/NavierStokesSpecialBoundaryConditions.cpp")
    symlinkHAMeRSFile(path_NavierStokesErrorStatistics_cpp, path_NavierStokesErrorStatistics_cpp_default, "src/apps/Navier-Stokes/NavierStokesErrorStatistics.cpp")

    # Link files related to flow models.
    symlinkHAMeRSFile(path_FlowModelSpecialSourceTerms_cpp, path_FlowModelSpecialSourceTerms_cpp_default, "src/flow/flow_models/FlowModelSpecialSourceTerms.cpp")
    symlinkHAMeRSFile(path_FlowModelStatisticsUtilitiesSingleSpecies_cpp, path_FlowModelStatisticsUtilitiesSingleSpecies_cpp_default, "src/flow/flow_models/single-species/FlowModelStatisticsUtilitiesSingleSpecies.cpp")
    symlinkHAMeRSFile(path_FlowModelStatisticsUtilitiesFourEqnConservative_cpp, path_FlowModelStatisticsUtilitiesFourEqnConservative_cpp_default, "src/flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.cpp")
    symlinkHAMeRSFile(path_FlowModelStatisticsUtilitiesFiveEqnAllaire_cpp, path_FlowModelStatisticsUtilitiesFiveEqnAllaire_cpp_default, "src/flow/flow_models/five-eqn_Allaire/FlowModelStatisticsUtilitiesFiveEqnAllaire.cpp")

    # Link files related to immersed boundaries.
    symlinkHAMeRSFile(path_ImmersedBoundaries_cpp, path_ImmersedBoundaries_cpp_default, "src/util/immersed_boundaries/ImmersedBoundaries.cpp")

