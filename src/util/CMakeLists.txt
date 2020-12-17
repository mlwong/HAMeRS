# Define a variable utils_source_files containing
# a list of the source files for the utils library
set(utils_source_files
    empty.cpp
)

# Create a library called utils which includes the 
# source files defined in utils_source_files
add_library(utils ${utils_source_files})

TARGET_LINK_LIBRARIES(utils derivatives differences filters
  gradient_sensors wavelet_transform mixing_rules
  basic_boundary_conditions MPI_helpers)
