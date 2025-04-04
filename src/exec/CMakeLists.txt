# Compile main from main.cpp
add_executable(main main.cpp)
add_executable(sphere sphere.cpp)

# main should link to the Euler and Navier-Stokes library
target_link_libraries(main Euler Navier_Stokes visit_data_writer gfortran)
target_link_libraries(sphere tecio)
