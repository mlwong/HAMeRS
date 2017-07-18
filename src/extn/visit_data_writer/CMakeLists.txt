# Define a visit_data_writer_source_files containing
# a list of the source files for the visit_data_writer library
set(visit_data_writer_source_files
    ExtendedVisItDataWriter.cpp
)

# Create a library called visit_data_writer which includes the 
# source files defined in visit_data_writer_source_files
add_library(visit_data_writer ${visit_data_writer_source_files})

TARGET_LINK_LIBRARIES(visit_data_writer SAMRAI_appu SAMRAI_geom SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
