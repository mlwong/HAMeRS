# Define a variable MPI_helpers_source_files containing
# a list of the source files for the MPI helpers library
set(MPI_helpers_source_files
    MPIHelper.cpp
    MPIHelperMaxMin.cpp
    MPIHelperAverage.cpp
    MPIHelperCorrelation.cpp
    MPIHelperCentroid.cpp
    )

# Create a library called MPI_helpers which includes the 
# source files defined in MPI_helpers_source_files
add_library(MPI_helpers ${MPI_helpers_source_files})

TARGET_LINK_LIBRARIES(MPI_helpers SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
