# Define a variable feature_driven_tagger_source_files containing
# a list of the source files for the feature_driven_tagger library
set(feature_driven_tagger_source_files
    FeatureDrivenTagger.cpp
)

# Create a library called feature_driven_tagger which includes the 
# source files defined in feature_driven_tagger_source_files
add_library(feature_driven_tagger ${feature_driven_tagger_source_files})

TARGET_LINK_LIBRARIES(feature_driven_tagger
  equation_of_state algs_fortran
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  hdf5 gfortran)
