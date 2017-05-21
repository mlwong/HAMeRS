# Define a patch_hierarchies_source_files containing
# a list of the source files for the patch_hierarchies library
set(patch_hierarchies_source_files
    ExtendedUncoveredBoxIterator.cpp
    ExtendedFlattenedHierarchy.cpp
)

# Create a library called patch_hierarchies which includes the 
# source files defined in patch_hierarchies_source_files
add_library(patch_hierarchies ${patch_hierarchies_source_files})

TARGET_LINK_LIBRARIES(patch_hierarchies SAMRAI_hier SAMRAI_tbox)
