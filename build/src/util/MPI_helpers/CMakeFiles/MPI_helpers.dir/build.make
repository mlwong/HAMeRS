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
include src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/depend.make

# Include the progress variables for this target.
include src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/progress.make

# Include the compile flags for this target's objects.
include src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.o: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.o: ../src/util/MPI_helpers/MPIHelper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelper.cpp

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelper.cpp > CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.i

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelper.cpp -o CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.s

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.o: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.o: ../src/util/MPI_helpers/MPIHelperGrid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperGrid.cpp

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperGrid.cpp > CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.i

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperGrid.cpp -o CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.s

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.o: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.o: ../src/util/MPI_helpers/MPIHelperMaxMin.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperMaxMin.cpp

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperMaxMin.cpp > CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.i

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperMaxMin.cpp -o CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.s

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.o: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.o: ../src/util/MPI_helpers/MPIHelperAverage.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperAverage.cpp

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperAverage.cpp > CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.i

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperAverage.cpp -o CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.s

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.o: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.o: ../src/util/MPI_helpers/MPIHelperCorrelation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperCorrelation.cpp

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperCorrelation.cpp > CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.i

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperCorrelation.cpp -o CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.s

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.o: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/flags.make
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.o: ../src/util/MPI_helpers/MPIHelperCentroid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.o"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.o -c /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperCentroid.cpp

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.i"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperCentroid.cpp > CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.i

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.s"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && mpiicpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers/MPIHelperCentroid.cpp -o CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.s

# Object files for target MPI_helpers
MPI_helpers_OBJECTS = \
"CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.o" \
"CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.o" \
"CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.o" \
"CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.o" \
"CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.o" \
"CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.o"

# External object files for target MPI_helpers
MPI_helpers_EXTERNAL_OBJECTS =

src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelper.cpp.o
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperGrid.cpp.o
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperMaxMin.cpp.o
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperAverage.cpp.o
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCorrelation.cpp.o
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/MPIHelperCentroid.cpp.o
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/build.make
src/util/MPI_helpers/libMPI_helpers.a: src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/akula/HAMeRS_IB/HAMeRS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library libMPI_helpers.a"
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && $(CMAKE_COMMAND) -P CMakeFiles/MPI_helpers.dir/cmake_clean_target.cmake
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MPI_helpers.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/build: src/util/MPI_helpers/libMPI_helpers.a

.PHONY : src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/build

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/clean:
	cd /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers && $(CMAKE_COMMAND) -P CMakeFiles/MPI_helpers.dir/cmake_clean.cmake
.PHONY : src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/clean

src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/depend:
	cd /home/akula/HAMeRS_IB/HAMeRS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/akula/HAMeRS_IB/HAMeRS /home/akula/HAMeRS_IB/HAMeRS/src/util/MPI_helpers /home/akula/HAMeRS_IB/HAMeRS/build /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers /home/akula/HAMeRS_IB/HAMeRS/build/src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/util/MPI_helpers/CMakeFiles/MPI_helpers.dir/depend

