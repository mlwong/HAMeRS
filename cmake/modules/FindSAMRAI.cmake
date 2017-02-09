# Try to find SAMRAI library
# Once done this will define
#  SAMRAI_FOUND - if system found SAMRAI library
#  SAMRAI_INCLUDE_DIRS - The SAMRAI include directories
#  SAMRAI_LIBRARIES - The libraries needed to use SAMRAI
#  SAMRAI_DEFINITIONS - Compiler switches required for using SAMRAI

set (SAMRAI_ROOT_DIR $ENV{SAMRAI_ROOT})
# Uncomment the following line to print which directory CMake is looking in.
MESSAGE(STATUS "SAMRAI_ROOT_DIR: " ${SAMRAI_ROOT_DIR})

find_path(SAMRAI_INCLUDE_DIR
    NAMES SAMRAI/SAMRAI_config.h
    PATHS ${SAMRAI_ROOT_DIR}/include)

find_path(SAMRAI_LIBRARY_DIRS
    NAMES libSAMRAI_tbox.a
    PATHS ${SAMRAI_ROOT_DIR}/lib
)

if(SAMRAI_INCLUDE_DIR AND SAMRAI_LIBRARY_DIRS)
    set(SAMRAI_FOUND TRUE)
endif(SAMRAI_INCLUDE_DIR AND SAMRAI_LIBRARY_DIRS)

if (SAMRAI_FOUND)
    if (NOT SAMRAI_FIND_QUIETLY)
        message(STATUS "Found SAMRAI: ${SAMRAI_LIBRARY_DIRS}")
    endif (NOT SAMRAI_FIND_QUIETLY)
    else (SAMRAI_FOUND)
    if (SAMRAI_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find SAMRAI")
    endif (SAMRAI_FIND_REQUIRED)
endif (SAMRAI_FOUND)

