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
include src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/depend.make

# Include the progress variables for this target.
include src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/progress.make

# Include the compile flags for this target's objects.
include src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/flags.make

src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.o: src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/flags.make
src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.o: ../src/util/immersed_boundaries/ImmersedBoundaries.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/immersed_boundaries/ImmersedBoundaries.cpp

src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/immersed_boundaries/ImmersedBoundaries.cpp > CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.i

src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/immersed_boundaries/ImmersedBoundaries.cpp -o CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.s

# Object files for target immersed_boundaries
immersed_boundaries_OBJECTS = \
"CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.o"

# External object files for target immersed_boundaries
immersed_boundaries_EXTERNAL_OBJECTS =

src/util/immersed_boundaries/libimmersed_boundaries.a: src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/ImmersedBoundaries.cpp.o
src/util/immersed_boundaries/libimmersed_boundaries.a: src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/build.make
src/util/immersed_boundaries/libimmersed_boundaries.a: src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libimmersed_boundaries.a"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries && $(CMAKE_COMMAND) -P CMakeFiles/immersed_boundaries.dir/cmake_clean_target.cmake
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/immersed_boundaries.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/build: src/util/immersed_boundaries/libimmersed_boundaries.a

.PHONY : src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/build

src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/clean:
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries && $(CMAKE_COMMAND) -P CMakeFiles/immersed_boundaries.dir/cmake_clean.cmake
.PHONY : src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/clean

src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/depend:
	cd /home/akula/HAMeRS_IB/HAMeRS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/akula/HAMeRS_IB/HAMeRS /home/akula/HAMeRS_IB/HAMeRS/src/util/immersed_boundaries /home/akula/HAMeRS_IB/HAMeRS/build /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries /home/akula/HAMeRS_IB/HAMeRS/build/src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/util/immersed_boundaries/CMakeFiles/immersed_boundaries.dir/depend

