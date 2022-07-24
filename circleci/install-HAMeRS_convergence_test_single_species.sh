cd src/apps/Euler
ln -sf ${HAMERS_ROOT}/problems/Euler/initial_conditions/ConvergenceSingleSpecies.cpp EulerInitialConditions.cpp 
ln -sf ${HAMERS_ROOT}/problems/Euler/error_statistics/ConvergenceSingleSpecies.cpp EulerErrorStatistics.cpp
cd ../../..
mkdir build_convergence_test_single_species
cd build_convergence_test_single_species
export SAMRAI_ROOT=$SAMRAI_ROOT_NO_BOOST; cmake .. -DHAMERS_USE_BOOST=OFF -DHAMERS_ENABLE_SIMD=OFF
make -j 4
cd ..
