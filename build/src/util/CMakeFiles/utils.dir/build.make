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
include src/util/CMakeFiles/utils.dir/depend.make

# Include the progress variables for this target.
include src/util/CMakeFiles/utils.dir/progress.make

# Include the compile flags for this target's objects.
include src/util/CMakeFiles/utils.dir/flags.make

src/util/CMakeFiles/utils.dir/empty.cpp.o: src/util/CMakeFiles/utils.dir/flags.make
src/util/CMakeFiles/utils.dir/empty.cpp.o: ../src/util/empty.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/util/CMakeFiles/utils.dir/empty.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/empty.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/empty.cpp

src/util/CMakeFiles/utils.dir/empty.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/empty.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/empty.cpp > CMakeFiles/utils.dir/empty.cpp.i

src/util/CMakeFiles/utils.dir/empty.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/empty.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/empty.cpp -o CMakeFiles/utils.dir/empty.cpp.s

# Object files for target utils
utils_OBJECTS = \
"CMakeFiles/utils.dir/empty.cpp.o"

# External object files for target utils
utils_EXTERNAL_OBJECTS =

src/util/libutils.a: src/util/CMakeFiles/utils.dir/empty.cpp.o
src/util/libutils.a: src/util/CMakeFiles/utils.dir/build.make
src/util/libutils.a: src/util/CMakeFiles/utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libutils.a"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean_target.cmake
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/util/CMakeFiles/utils.dir/build: src/util/libutils.a

.PHONY : src/util/CMakeFiles/utils.dir/build

src/util/CMakeFiles/utils.dir/clean:
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean.cmake
.PHONY : src/util/CMakeFiles/utils.dir/clean

src/util/CMakeFiles/utils.dir/depend:
	cd /home/akula/HAMeRS_IB/HAMeRS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/akula/HAMeRS_IB/HAMeRS /home/akula/HAMeRS_IB/HAMeRS/src/util /home/akula/HAMeRS_IB/HAMeRS/build /home/akula/HAMeRS_IB/HAMeRS/build/src/util /home/akula/HAMeRS_IB/HAMeRS/build/src/util/CMakeFiles/utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/util/CMakeFiles/utils.dir/depend

