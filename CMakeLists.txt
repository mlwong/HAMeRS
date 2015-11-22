cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER gfortran)

project(AMR-WCNS C CXX Fortran)

find_package(MPI)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

# Directories to look for include files
include_directories(${AMR-WCNS_SOURCE_DIR}/include
                    ${AMR-WCNS_SOURCE_DIR}/hdf5/include
                    ${AMR-WCNS_SOURCE_DIR}/SAMRAI/include)

# Directories to look for libraries
link_directories(${AMR-WCNS_SOURCE_DIR}/hdf5/lib
                 ${AMR-WCNS_SOURCE_DIR}/SAMRAI/lib
		 ${AMR-WCNS_BINARY_DIR}/patch_strategy
		 ${AMR-WCNS_BINARY_DIR}/integrator/fortran
                 ${AMR-WCNS_BINARY_DIR}/integrator
		 ${AMR-WCNS_BINARY_DIR}/applications/Euler
                 ${AMR-WCNS_BINARY_DIR}/flow_model/boundary_conditions/Euler
		 ${AMR-WCNS_BINARY_DIR}/flow_model/convective_flux_reconstructor
		 ${AMR-WCNS_BINARY_DIR}/flow_model/equation_of_state
		 ${AMR-WCNS_BINARY_DIR}/flow_model/feature_driven_tagger
                 ${AMR-WCNS_BINARY_DIR}/flow_model/initial_conditions
		 ${AMR-WCNS_BINARY_DIR}/flow_model/Riemann_solver)

add_subdirectory(src/integrator/fortran)

# Add the compile flags
add_definitions("-std=c++11 -g -O3 -Wall -Wno-deprecated")

add_subdirectory(src/patch_strategy)
add_subdirectory(src/integrator)
add_subdirectory(src/applications/Euler)
add_subdirectory(src/exec)
add_subdirectory(src/flow_model/boundary_conditions/Euler)
add_subdirectory(src/flow_model/convective_flux_reconstructor)
add_subdirectory(src/flow_model/equation_of_state)
add_subdirectory(src/flow_model/feature_driven_tagger)
add_subdirectory(src/flow_model/initial_conditions)
add_subdirectory(src/flow_model/Riemann_solver)
