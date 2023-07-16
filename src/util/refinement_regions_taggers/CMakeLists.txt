# Define a variable refine_regions_taggers_source_files containing
# a list of the source files for the refine_regions_taggers library
set(refine_regions_taggers_source_files
    RefineRegionsTagger.cpp
    )

# Create a library called refine_regions_taggers which includes the 
# source files defined in refine_regions_taggers_source_files
add_library(refine_regions_taggers ${refine_regions_taggers_source_files})

TARGET_LINK_LIBRARIES(refine_regions_taggers SAMRAI_geom SAMRAI_hier SAMRAI_tbox)
