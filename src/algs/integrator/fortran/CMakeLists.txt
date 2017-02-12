# Define a variable application_source_files containing a list of the source
# files for the algs_fortran library
set(algs_fortran_source_files
    algs_upfluxsum3d.m4
    algs_upfluxsum3d.f
    algs_upfluxsum2d.m4
    algs_upfluxsum2d.f
    algs_upfluxsum1d.m4
    algs_upfluxsum1d.f
    algs_nodeouternodeops3d.m4
    algs_nodeouternodeops3d.f
    algs_nodeouternodeops2d.m4
    algs_nodeouternodeops2d.f
)

# Create a library called algs_fortran which includes the source files defined
# in algs_fortran_source_files
add_library(algs_fortran ${algs_fortran_source_files})

