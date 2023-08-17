# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /share/apps/cmake/3.12.4/bin/cmake

# The command to remove a file.
RM = /share/apps/cmake/3.12.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/akula/HAMeRS_IB/HAMeRS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/akula/HAMeRS_IB/HAMeRS/build

# Include any dependencies generated for this target.
include src/util/differences/CMakeFiles/differences.dir/depend.make

# Include the progress variables for this target.
include src/util/differences/CMakeFiles/differences.dir/progress.make

# Include the compile flags for this target's objects.
include src/util/differences/CMakeFiles/differences.dir/flags.make

src/util/differences/CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.o: src/util/differences/CMakeFiles/differences.dir/flags.make
src/util/differences/CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.o: ../src/util/differences/DifferenceFirstOrder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/util/differences/CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/differences/DifferenceFirstOrder.cpp

src/util/differences/CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/differences/DifferenceFirstOrder.cpp > CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.i

src/util/differences/CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/differences/DifferenceFirstOrder.cpp -o CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.s

src/util/differences/CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.o: src/util/differences/CMakeFiles/differences.dir/flags.make
src/util/differences/CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.o: ../src/util/differences/DifferenceSecondOrder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/util/differences/CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/differences/DifferenceSecondOrder.cpp

src/util/differences/CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/differences/DifferenceSecondOrder.cpp > CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.i

src/util/differences/CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/differences/DifferenceSecondOrder.cpp -o CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.s

# Object files for target differences
differences_OBJECTS = \
"CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.o" \
"CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.o"

# External object files for target differences
differences_EXTERNAL_OBJECTS =

src/util/differences/libdifferences.a: src/util/differences/CMakeFiles/differences.dir/DifferenceFirstOrder.cpp.o
src/util/differences/libdifferences.a: src/util/differences/CMakeFiles/differences.dir/DifferenceSecondOrder.cpp.o
src/util/differences/libdifferences.a: src/util/differences/CMakeFiles/differences.dir/build.make
src/util/differences/libdifferences.a: src/util/differences/CMakeFiles/differences.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libdifferences.a"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && $(CMAKE_COMMAND) -P CMakeFiles/differences.dir/cmake_clean_target.cmake
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/differences.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/util/differences/CMakeFiles/differences.dir/build: src/util/differences/libdifferences.a

.PHONY : src/util/differences/CMakeFiles/differences.dir/build

src/util/differences/CMakeFiles/differences.dir/clean:
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences && $(CMAKE_COMMAND) -P CMakeFiles/differences.dir/cmake_clean.cmake
.PHONY : src/util/differences/CMakeFiles/differences.dir/clean

src/util/differences/CMakeFiles/differences.dir/depend:
	cd /home/akula/HAMeRS_IB/HAMeRS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/akula/HAMeRS_IB/HAMeRS /home/akula/HAMeRS_IB/HAMeRS/src/util/differences /home/akula/HAMeRS_IB/HAMeRS/build /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences /home/akula/HAMeRS_IB/HAMeRS/build/src/util/differences/CMakeFiles/differences.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/util/differences/CMakeFiles/differences.dir/depend

