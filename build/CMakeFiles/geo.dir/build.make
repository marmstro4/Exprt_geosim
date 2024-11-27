# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/michaela/straw/MichaelGarfield/expert_geosim

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/michaela/straw/MichaelGarfield/expert_geosim/build

# Include any dependencies generated for this target.
include CMakeFiles/geo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/geo.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/geo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/geo.dir/flags.make

CMakeFiles/geo.dir/geo.cc.o: CMakeFiles/geo.dir/flags.make
CMakeFiles/geo.dir/geo.cc.o: /home/michaela/straw/MichaelGarfield/expert_geosim/geo.cc
CMakeFiles/geo.dir/geo.cc.o: CMakeFiles/geo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/michaela/straw/MichaelGarfield/expert_geosim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/geo.dir/geo.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/geo.dir/geo.cc.o -MF CMakeFiles/geo.dir/geo.cc.o.d -o CMakeFiles/geo.dir/geo.cc.o -c /home/michaela/straw/MichaelGarfield/expert_geosim/geo.cc

CMakeFiles/geo.dir/geo.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/geo.dir/geo.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/michaela/straw/MichaelGarfield/expert_geosim/geo.cc > CMakeFiles/geo.dir/geo.cc.i

CMakeFiles/geo.dir/geo.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/geo.dir/geo.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/michaela/straw/MichaelGarfield/expert_geosim/geo.cc -o CMakeFiles/geo.dir/geo.cc.s

# Object files for target geo
geo_OBJECTS = \
"CMakeFiles/geo.dir/geo.cc.o"

# External object files for target geo
geo_EXTERNAL_OBJECTS =

geo: CMakeFiles/geo.dir/geo.cc.o
geo: CMakeFiles/geo.dir/build.make
geo: /home/michaela/Downloads/garfield/garfieldpp/install/lib/libGarfield.so.0.3.0
geo: /home/michaela/Downloads/garfield/garfieldpp/install/lib/libHeed.so
geo: /home/michaela/Downloads/garfield/garfieldpp/install/lib/libGarfieldRandom.so
geo: /home/michaela/Downloads/root_build/lib/libGdml.so
geo: /home/michaela/Downloads/root_build/lib/libGeom.so
geo: /home/michaela/Downloads/root_build/lib/libXMLIO.so
geo: /home/michaela/Downloads/root_build/lib/libGraf3d.so
geo: /home/michaela/Downloads/root_build/lib/libGpad.so
geo: /home/michaela/Downloads/root_build/lib/libGraf.so
geo: /home/michaela/Downloads/root_build/lib/libHist.so
geo: /home/michaela/Downloads/root_build/lib/libMatrix.so
geo: /home/michaela/Downloads/root_build/lib/libMathCore.so
geo: /home/michaela/Downloads/root_build/lib/libImt.so
geo: /home/michaela/Downloads/root_build/lib/libMultiProc.so
geo: /home/michaela/Downloads/root_build/lib/libNet.so
geo: /home/michaela/Downloads/root_build/lib/libRIO.so
geo: /home/michaela/Downloads/root_build/lib/libThread.so
geo: /home/michaela/Downloads/root_build/lib/libCore.so
geo: /usr/lib/x86_64-linux-gnu/libgsl.so
geo: /usr/lib/x86_64-linux-gnu/libgslcblas.so
geo: /home/michaela/Downloads/garfield/garfieldpp/install/lib/libmagboltz.so.11
geo: /home/michaela/Downloads/garfield/garfieldpp/install/lib/libdegrade.so.3
geo: CMakeFiles/geo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/michaela/straw/MichaelGarfield/expert_geosim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable geo"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/geo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/geo.dir/build: geo
.PHONY : CMakeFiles/geo.dir/build

CMakeFiles/geo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/geo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/geo.dir/clean

CMakeFiles/geo.dir/depend:
	cd /home/michaela/straw/MichaelGarfield/expert_geosim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/michaela/straw/MichaelGarfield/expert_geosim /home/michaela/straw/MichaelGarfield/expert_geosim /home/michaela/straw/MichaelGarfield/expert_geosim/build /home/michaela/straw/MichaelGarfield/expert_geosim/build /home/michaela/straw/MichaelGarfield/expert_geosim/build/CMakeFiles/geo.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/geo.dir/depend

