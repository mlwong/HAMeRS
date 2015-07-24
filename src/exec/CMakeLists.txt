# Compile main from main.cpp
add_executable(main main.cpp)

# main should link to the Euler library
target_link_libraries(main Euler gfortran)
