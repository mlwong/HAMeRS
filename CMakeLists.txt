cmake_minimum_required(VERSION 2.8)

project(HAMeRS C CXX Fortran)

find_package(MPI)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

# Find the boost package
find_package(Boost 1.45 REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# Find the hdf5 package
find_package(HDF5 1.8.4 REQUIRED)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})

# Set the preprocessor variables
option(HAMERS_ENABLE_SIMD "HAMERS_ENABLE_SIMD" OFF)
option(HAMERS_DEBUG_CHECK_ASSERTIONS "HAMERS_DEBUG_CHECK_ASSERTIONS" ON)
option(HAMERS_DEBUG_CHECK_DEV_ASSERTIONS "HAMERS_DEBUG_CHECK_DEV_ASSERTIONS" ON)

configure_file(${HAMeRS_SOURCE_DIR}/include/HAMeRS_config.hpp.in ${HAMeRS_SOURCE_DIR}/include/HAMeRS_config.hpp)

# Set the compilers
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER gfortran)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Intel compiler suite
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3 -Wall -Wno-deprecated -xhost")
    elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3 -Wall -Wno-deprecated")
    endif()

# Standard GNU compilers
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3 -Wall -Wno-deprecated")
    elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -Wall -Wno-deprecated")
    endif()
endif()

# Directories to look for include files
include_directories(${HAMeRS_SOURCE_DIR}/include
                    /afs/ir.stanford.edu/users/w/o/wongml/AMR-WCNS/SAMRAI/include)

# Directories to look for libraries
link_directories(/afs/ir.stanford.edu/users/w/o/wongml/AMR-WCNS/SAMRAI/lib
                 ${HAMeRS_BINARY_DIR}/algs/patch_strategy
                 ${HAMeRS_BINARY_DIR}/algs/integrator/fortran
                 ${HAMeRS_BINARY_DIR}/algs/integrator
                 ${HAMeRS_BINARY_DIR}/apps/Euler
                 ${HAMeRS_BINARY_DIR}/apps/Navier-Stokes
                 ${HAMeRS_BINARY_DIR}/flow/flow_models
                 ${HAMeRS_BINARY_DIR}/flow/convective_flux_reconstructors
                 ${HAMeRS_BINARY_DIR}/flow/diffusive_flux_reconstructors
                 ${HAMeRS_BINARY_DIR}/flow/refinement_taggers
                 ${HAMeRS_BINARY_DIR}/util
                 ${HAMeRS_BINARY_DIR}/util/basic_boundary_conditions
                 ${HAMeRS_BINARY_DIR}/util/mixing_rules
                 ${HAMeRS_BINARY_DIR}/util/mixing_rules/equations_of_state)

add_subdirectory(src/algs/integrator/fortran)
add_subdirectory(src/algs/patch_strategy)
add_subdirectory(src/algs/integrator)
add_subdirectory(src/apps/Euler)
add_subdirectory(src/apps/Navier-Stokes)
add_subdirectory(src/exec)
add_subdirectory(src/flow/flow_models)
add_subdirectory(src/flow/convective_flux_reconstructors)
add_subdirectory(src/flow/diffusive_flux_reconstructors)
add_subdirectory(src/flow/refinement_taggers)
add_subdirectory(src/util)
add_subdirectory(src/util/basic_boundary_conditions)
add_subdirectory(src/util/mixing_rules)
add_subdirectory(src/util/mixing_rules/equations_of_state)
