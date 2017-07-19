# Define a variable refinement_taggers_source_files containing
# a list of the source files for the refinement_taggers library
set(refinement_taggers_source_files
    ValueTagger.cpp
    GradientTagger.cpp
    MultiresolutionTagger.cpp
)

# Create a library called refinement_taggers which includes the 
# source files defined in refinement_taggers_source_files
add_library(refinement_taggers ${refinement_taggers_source_files})

TARGET_LINK_LIBRARIES(refinement_taggers utils flow_models
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
