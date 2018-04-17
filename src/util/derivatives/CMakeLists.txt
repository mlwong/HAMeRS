# Define a variable derivatives_source_files containing
# a list of the source files for the derivatives library
set(derivatives_source_files
    Derivative.cpp
    DerivativeFirstOrder.cpp
    DerivativeSecondOrder.cpp
    )

# Create a library called derivatives which includes the 
# source files defined in derivatives_source_files
add_library(derivatives ${derivatives_source_files})

TARGET_LINK_LIBRARIES(derivatives SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
