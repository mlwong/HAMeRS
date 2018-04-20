# Define a variable wavelet_transform_source_files containing
# a list of the source files for the wavelet_transform library
set(wavelet_transform_source_files
    WaveletTransformHarten.cpp
    )

# Create a library called wavelet_transform which includes the 
# source files defined in wavelet_transform_source_files
add_library(wavelet_transform ${wavelet_transform_source_files})

TARGET_LINK_LIBRARIES(wavelet_transform SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
