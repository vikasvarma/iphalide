# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/vikasvarma/Documents/Development/iphalide

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/vikasvarma/Documents/Development/iphalide/build

# Include any dependencies generated for this target.
include CMakeFiles/conv_halide.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/conv_halide.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/conv_halide.dir/flags.make

CMakeFiles/conv_halide.dir/src/conv.cpp.o: CMakeFiles/conv_halide.dir/flags.make
CMakeFiles/conv_halide.dir/src/conv.cpp.o: ../src/conv.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/vikasvarma/Documents/Development/iphalide/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/conv_halide.dir/src/conv.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/conv_halide.dir/src/conv.cpp.o -c /Users/vikasvarma/Documents/Development/iphalide/src/conv.cpp

CMakeFiles/conv_halide.dir/src/conv.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/conv_halide.dir/src/conv.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/vikasvarma/Documents/Development/iphalide/src/conv.cpp > CMakeFiles/conv_halide.dir/src/conv.cpp.i

CMakeFiles/conv_halide.dir/src/conv.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/conv_halide.dir/src/conv.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/vikasvarma/Documents/Development/iphalide/src/conv.cpp -o CMakeFiles/conv_halide.dir/src/conv.cpp.s

CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.o: CMakeFiles/conv_halide.dir/flags.make
CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.o: /usr/local/share/tools/GenGen.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/vikasvarma/Documents/Development/iphalide/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.o -c /usr/local/share/tools/GenGen.cpp

CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /usr/local/share/tools/GenGen.cpp > CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.i

CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /usr/local/share/tools/GenGen.cpp -o CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.s

# Object files for target conv_halide
conv_halide_OBJECTS = \
"CMakeFiles/conv_halide.dir/src/conv.cpp.o" \
"CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.o"

# External object files for target conv_halide
conv_halide_EXTERNAL_OBJECTS =

conv_halide: CMakeFiles/conv_halide.dir/src/conv.cpp.o
conv_halide: CMakeFiles/conv_halide.dir/usr/local/share/tools/GenGen.cpp.o
conv_halide: CMakeFiles/conv_halide.dir/build.make
conv_halide: /usr/local/lib/libHalide.12.0.1.dylib
conv_halide: CMakeFiles/conv_halide.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/vikasvarma/Documents/Development/iphalide/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable conv_halide"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/conv_halide.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/conv_halide.dir/build: conv_halide

.PHONY : CMakeFiles/conv_halide.dir/build

CMakeFiles/conv_halide.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/conv_halide.dir/cmake_clean.cmake
.PHONY : CMakeFiles/conv_halide.dir/clean

CMakeFiles/conv_halide.dir/depend:
	cd /Users/vikasvarma/Documents/Development/iphalide/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/vikasvarma/Documents/Development/iphalide /Users/vikasvarma/Documents/Development/iphalide /Users/vikasvarma/Documents/Development/iphalide/build /Users/vikasvarma/Documents/Development/iphalide/build /Users/vikasvarma/Documents/Development/iphalide/build/CMakeFiles/conv_halide.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/conv_halide.dir/depend

