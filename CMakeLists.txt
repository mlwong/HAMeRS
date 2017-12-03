cmake_minimum_required(VERSION 2.8)

project(HAMeRS C CXX Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)


find_package(MPI)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Find the boost package
find_package(Boost 1.45 REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

# Find the hdf5 package
find_package(HDF5 1.8.4 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})
link_directories(${HDF5_LIBRARY_DIRS})

# Find the zlib package
if (WITH_ZLIB)
    set(ZLIB_LIBRARIES $ENV{ZLIB_ROOT})
    set(ZLIB_LIBRARY_DIRS "${ZLIB_ROOT}/lib")
    set(ZLIB_INCLUDE_DIRS "${ZLIB_ROOT}/include")
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARY_DIRS})
endif()

# Find the SAMRAI package
find_package(SAMRAI REQUIRED)
include_directories(${SAMRAI_INCLUDE_DIR})
link_directories(${SAMRAI_LIBRARY_DIRS})

message(STATUS "SAMRAI_INCLUDE_DIR: " ${SAMRAI_INCLUDE_DIR})

# Set the preprocessor variables
option(HAMERS_ENABLE_SIMD "HAMERS_ENABLE_SIMD" OFF)
option(HAMERS_DEBUG_CHECK_ASSERTIONS "HAMERS_DEBUG_CHECK_ASSERTIONS" ON)
option(HAMERS_DEBUG_CHECK_DEV_ASSERTIONS "HAMERS_DEBUG_CHECK_DEV_ASSERTIONS" ON)

configure_file(${HAMeRS_SOURCE_DIR}/include/HAMeRS_config.hpp.in ${HAMeRS_SOURCE_DIR}/include/HAMeRS_config.hpp)

# Set the compilers
set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})
set(CMAKE_Fortran_COMPILER $ENV{F77})

message(STATUS "CMAKE_C_COMPILER_ID: " ${CMAKE_C_COMPILER_ID})
message(STATUS "CMAKE_CXX_COMPILER_ID: " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "CMAKE_Fortran_COMPILER_ID: " ${CMAKE_Fortran_COMPILER_ID})

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Intel compiler suite
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wno-deprecated -qopenmp -qoverride-limits -xhost -qopt-report")
    elseif (CMAKE_BUILD_TYPE MATCHES "KNL")
        set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wno-deprecated -qopenmp -axMIC-AVX512 -qopt-report")
    elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -Wall -Wno-deprecated")
    endif()

# Standard GNU compilers
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wno-deprecated -fopenmp")
    elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -Wall -Wno-deprecated")
    endif()

# Clang compilers
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wno-deprecated -fopenmp")
    elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -Wall -Wno-deprecated")
    endif()

endif()

# Directories to look for include files
include_directories(${HAMeRS_SOURCE_DIR}/include)

# Directories to look for libraries
link_directories(${HAMeRS_BINARY_DIR}/algs/patch_strategy
                 ${HAMeRS_BINARY_DIR}/algs/integrator
                 ${HAMeRS_BINARY_DIR}/apps/Euler
                 ${HAMeRS_BINARY_DIR}/apps/Navier-Stokes
                 ${HAMeRS_BINARY_DIR}/flow/flow_models
                 ${HAMeRS_BINARY_DIR}/flow/convective_flux_reconstructors
                 ${HAMeRS_BINARY_DIR}/flow/diffusive_flux_reconstructors
                 ${HAMeRS_BINARY_DIR}/flow/nonconservative_diffusive_flux_divergence_operators
                 ${HAMeRS_BINARY_DIR}/flow/refinement_taggers
                 ${HAMeRS_BINARY_DIR}/extn/patch_hierarchies
                 ${HAMeRS_BINARY_DIR}/extn/visit_data_writer
                 ${HAMeRS_BINARY_DIR}/util
                 ${HAMeRS_BINARY_DIR}/util/basic_boundary_conditions
                 ${HAMeRS_BINARY_DIR}/util/mixing_rules
                 ${HAMeRS_BINARY_DIR}/util/mixing_rules/equations_of_state)

add_subdirectory(src/algs/patch_strategy)
add_subdirectory(src/algs/integrator)
add_subdirectory(src/apps/Euler)
add_subdirectory(src/apps/Navier-Stokes)
add_subdirectory(src/exec)
add_subdirectory(src/flow/flow_models)
add_subdirectory(src/flow/convective_flux_reconstructors)
add_subdirectory(src/flow/diffusive_flux_reconstructors)
add_subdirectory(src/flow/nonconservative_diffusive_flux_divergence_operators)
add_subdirectory(src/flow/refinement_taggers)
add_subdirectory(src/extn/patch_hierarchies)
add_subdirectory(src/extn/visit_data_writer)
add_subdirectory(src/util)
add_subdirectory(src/util/basic_boundary_conditions)
add_subdirectory(src/util/mixing_rules)
add_subdirectory(src/util/mixing_rules/equations_of_state)
