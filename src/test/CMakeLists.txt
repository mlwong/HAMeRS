# Compile all tests
add_executable(test_mixing_rules test_mixing_rules.cpp)

# All tests should link to the Euler and Navier-Stokes library
target_link_libraries(test_mixing_rules Euler Navier_Stokes gfortran)
