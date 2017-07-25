# Define a variable gradient_sensors_source_files containing
# a list of the source files for the gradient_sensors library
set(gradient_sensors_source_files
    GradientSensorJameson.cpp
    )

# Create a library called gradient_sensors which includes the 
# source files defined in gradient_sensors_source_files
add_library(gradient_sensors ${gradient_sensors_source_files})

TARGET_LINK_LIBRARIES(gradient_sensors SAMRAI_geom SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
