# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /cygdrive/c/Users/Darren/.CLion2018.1/system/cygwin_cmake/bin/cmake.exe

# The command to remove a file.
RM = /cygdrive/c/Users/Darren/.CLion2018.1/system/cygwin_cmake/bin/cmake.exe -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cygdrive/c/Git/MinhashSketch

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cygdrive/c/Git/MinhashSketch/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MinhashSketch.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MinhashSketch.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MinhashSketch.dir/flags.make

CMakeFiles/MinhashSketch.dir/main.cpp.o: CMakeFiles/MinhashSketch.dir/flags.make
CMakeFiles/MinhashSketch.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Git/MinhashSketch/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MinhashSketch.dir/main.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MinhashSketch.dir/main.cpp.o -c /cygdrive/c/Git/MinhashSketch/main.cpp

CMakeFiles/MinhashSketch.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MinhashSketch.dir/main.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Git/MinhashSketch/main.cpp > CMakeFiles/MinhashSketch.dir/main.cpp.i

CMakeFiles/MinhashSketch.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MinhashSketch.dir/main.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Git/MinhashSketch/main.cpp -o CMakeFiles/MinhashSketch.dir/main.cpp.s

CMakeFiles/MinhashSketch.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/MinhashSketch.dir/main.cpp.o.requires

CMakeFiles/MinhashSketch.dir/main.cpp.o.provides: CMakeFiles/MinhashSketch.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/MinhashSketch.dir/build.make CMakeFiles/MinhashSketch.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/MinhashSketch.dir/main.cpp.o.provides

CMakeFiles/MinhashSketch.dir/main.cpp.o.provides.build: CMakeFiles/MinhashSketch.dir/main.cpp.o


CMakeFiles/MinhashSketch.dir/Minhash.cpp.o: CMakeFiles/MinhashSketch.dir/flags.make
CMakeFiles/MinhashSketch.dir/Minhash.cpp.o: ../Minhash.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Git/MinhashSketch/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MinhashSketch.dir/Minhash.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MinhashSketch.dir/Minhash.cpp.o -c /cygdrive/c/Git/MinhashSketch/Minhash.cpp

CMakeFiles/MinhashSketch.dir/Minhash.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MinhashSketch.dir/Minhash.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Git/MinhashSketch/Minhash.cpp > CMakeFiles/MinhashSketch.dir/Minhash.cpp.i

CMakeFiles/MinhashSketch.dir/Minhash.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MinhashSketch.dir/Minhash.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Git/MinhashSketch/Minhash.cpp -o CMakeFiles/MinhashSketch.dir/Minhash.cpp.s

CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.requires:

.PHONY : CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.requires

CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.provides: CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.requires
	$(MAKE) -f CMakeFiles/MinhashSketch.dir/build.make CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.provides.build
.PHONY : CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.provides

CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.provides.build: CMakeFiles/MinhashSketch.dir/Minhash.cpp.o


CMakeFiles/MinhashSketch.dir/Utils.cpp.o: CMakeFiles/MinhashSketch.dir/flags.make
CMakeFiles/MinhashSketch.dir/Utils.cpp.o: ../Utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Git/MinhashSketch/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MinhashSketch.dir/Utils.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MinhashSketch.dir/Utils.cpp.o -c /cygdrive/c/Git/MinhashSketch/Utils.cpp

CMakeFiles/MinhashSketch.dir/Utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MinhashSketch.dir/Utils.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Git/MinhashSketch/Utils.cpp > CMakeFiles/MinhashSketch.dir/Utils.cpp.i

CMakeFiles/MinhashSketch.dir/Utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MinhashSketch.dir/Utils.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Git/MinhashSketch/Utils.cpp -o CMakeFiles/MinhashSketch.dir/Utils.cpp.s

CMakeFiles/MinhashSketch.dir/Utils.cpp.o.requires:

.PHONY : CMakeFiles/MinhashSketch.dir/Utils.cpp.o.requires

CMakeFiles/MinhashSketch.dir/Utils.cpp.o.provides: CMakeFiles/MinhashSketch.dir/Utils.cpp.o.requires
	$(MAKE) -f CMakeFiles/MinhashSketch.dir/build.make CMakeFiles/MinhashSketch.dir/Utils.cpp.o.provides.build
.PHONY : CMakeFiles/MinhashSketch.dir/Utils.cpp.o.provides

CMakeFiles/MinhashSketch.dir/Utils.cpp.o.provides.build: CMakeFiles/MinhashSketch.dir/Utils.cpp.o


# Object files for target MinhashSketch
MinhashSketch_OBJECTS = \
"CMakeFiles/MinhashSketch.dir/main.cpp.o" \
"CMakeFiles/MinhashSketch.dir/Minhash.cpp.o" \
"CMakeFiles/MinhashSketch.dir/Utils.cpp.o"

# External object files for target MinhashSketch
MinhashSketch_EXTERNAL_OBJECTS =

MinhashSketch.exe: CMakeFiles/MinhashSketch.dir/main.cpp.o
MinhashSketch.exe: CMakeFiles/MinhashSketch.dir/Minhash.cpp.o
MinhashSketch.exe: CMakeFiles/MinhashSketch.dir/Utils.cpp.o
MinhashSketch.exe: CMakeFiles/MinhashSketch.dir/build.make
MinhashSketch.exe: CMakeFiles/MinhashSketch.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cygdrive/c/Git/MinhashSketch/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable MinhashSketch.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MinhashSketch.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MinhashSketch.dir/build: MinhashSketch.exe

.PHONY : CMakeFiles/MinhashSketch.dir/build

CMakeFiles/MinhashSketch.dir/requires: CMakeFiles/MinhashSketch.dir/main.cpp.o.requires
CMakeFiles/MinhashSketch.dir/requires: CMakeFiles/MinhashSketch.dir/Minhash.cpp.o.requires
CMakeFiles/MinhashSketch.dir/requires: CMakeFiles/MinhashSketch.dir/Utils.cpp.o.requires

.PHONY : CMakeFiles/MinhashSketch.dir/requires

CMakeFiles/MinhashSketch.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MinhashSketch.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MinhashSketch.dir/clean

CMakeFiles/MinhashSketch.dir/depend:
	cd /cygdrive/c/Git/MinhashSketch/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cygdrive/c/Git/MinhashSketch /cygdrive/c/Git/MinhashSketch /cygdrive/c/Git/MinhashSketch/cmake-build-debug /cygdrive/c/Git/MinhashSketch/cmake-build-debug /cygdrive/c/Git/MinhashSketch/cmake-build-debug/CMakeFiles/MinhashSketch.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MinhashSketch.dir/depend

