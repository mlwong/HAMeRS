# Install script for directory: /home/akula/HAMeRS_IB/HAMeRS

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "AMD")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/algs/patch_strategy/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/algs/integrator/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/apps/Euler/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/apps/Navier-Stokes/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/exec/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/flow/flow_models/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/flow/convective_flux_reconstructors/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/flow/diffusive_flux_reconstructors/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/flow/nonconservative_diffusive_flux_divergence_operators/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/flow/refinement_taggers/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/extn/patch_hierarchies/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/extn/visit_data_writer/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/test/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/basic_boundary_conditions/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/derivatives/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/filters/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/gradient_sensors/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/mixing_rules/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/mixing_rules/equations_of_state/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/refinement_regions_taggers/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/wavelet_transform/cmake_install.cmake")
  include("/home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/akula/HAMeRS_IB/HAMeRS/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
