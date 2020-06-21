# Define a variable filters_source_files containing
# a list of the source files for the filters library
set(filters_source_files
    Filter.cpp
    FilterNone.cpp
    FilterTruncatedGaussian.cpp
    )

# Create a library called filters which includes the 
# source files defined in filters_source_files
add_library(filters ${filters_source_files})

TARGET_LINK_LIBRARIES(filters SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
