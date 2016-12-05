# Compile main from main.cpp
add_executable(main main.cpp)

# main should link to the Euler and Navier-Stokes library
target_link_libraries(main Euler Navier_Stokes gfortran)

# Compile test from test.cpp
add_executable(test test.cpp)

# test should link to the Euler and Navier-Stokes library
target_link_libraries(test Euler Navier_Stokes gfortran)