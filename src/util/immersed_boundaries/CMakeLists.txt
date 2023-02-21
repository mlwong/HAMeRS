# Define a variable immersed_boundaries_source_files containing
# a list of the source files for the immersed_boundaries library
set(immersed_boundaries_source_files
    ImmersedBoundaries.cpp
    )

# Create a library called immersed_boundaries which includes the 
# source files defined in immersed_boundaries_source_files
add_library(immersed_boundaries ${immersed_boundaries_source_files})

TARGET_LINK_LIBRARIES(immersed_boundaries SAMRAI_geom SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
