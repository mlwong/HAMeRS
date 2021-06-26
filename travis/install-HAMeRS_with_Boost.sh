mkdir build_with_boost
cd build_with_boost
export SAMRAI_ROOT=$SAMRAI_ROOT_WITH_BOOST; cmake .. -DHAMERS_USE_BOOST=ON -DHAMERS_ENABLE_SIMD=OFF
make
cd ..
