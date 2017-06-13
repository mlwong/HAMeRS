# Define a transfer_source_files containing
# a list of the source files for the transfer library
set(transfer_source_files
    ExtendedRefineSchedule.cpp
    ExtendedRefineAlgorithm.cpp
)

# Create a library called transfer which includes the 
# source files defined in transfer_source_files
add_library(transfer ${transfer_source_files})

TARGET_LINK_LIBRARIES(transfer SAMRAI_xfer SAMRAI_hier SAMRAI_tbox)
