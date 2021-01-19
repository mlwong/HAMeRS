mkdir build_with_boost
cd build_with_boost
export SAMRAI_ROOT=$SAMRAI_ROOT_WITH_BOOST; cmake .. -DHAMERS_USE_BOOST=ON
make
cd ..