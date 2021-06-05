cd src/apps/Euler
ln -sf ${HAMERS_ROOT}/problems/Euler/initial_conditions/ConvergenceFiveEqnAllaire.cpp EulerInitialConditions.cpp 
ln -sf ${HAMERS_ROOT}/problems/Euler/error_statistics/ConvergenceFiveEqnAllaire.cpp EulerErrorStatistics.cpp
cd ../../..
mkdir build_convergence_test_five_eqn_allaire
cd build_convergence_test_five_eqn_allaire
export SAMRAI_ROOT=$SAMRAI_ROOT_NO_BOOST; cmake .. -DHAMERS_USE_BOOST=OFF -DHAMERS_ENABLE_SIMD=OFF
make
cd ..
