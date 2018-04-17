# Define a variable differences_source_files containing
# a list of the source files for the differences library
set(differences_source_files
    DifferenceFirstOrder.cpp
    DifferenceSecondOrder.cpp
    )

# Create a library called differences which includes the 
# source files defined in differences_source_files
add_library(differences ${differences_source_files})

TARGET_LINK_LIBRARIES(differences SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
