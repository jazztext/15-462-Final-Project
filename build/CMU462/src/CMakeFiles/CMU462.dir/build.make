# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.3

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zach/Documents/15462/asst3_pathtracer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zach/Documents/15462/asst3_pathtracer/build

# Include any dependencies generated for this target.
include CMU462/src/CMakeFiles/CMU462.dir/depend.make

# Include the progress variables for this target.
include CMU462/src/CMakeFiles/CMU462.dir/progress.make

# Include the compile flags for this target's objects.
include CMU462/src/CMakeFiles/CMU462.dir/flags.make

CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o: ../CMU462/src/vector2D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/vector2D.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector2D.cpp

CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/vector2D.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector2D.cpp > CMakeFiles/CMU462.dir/vector2D.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/vector2D.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector2D.cpp -o CMakeFiles/CMU462.dir/vector2D.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o: ../CMU462/src/vector3D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/vector3D.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector3D.cpp

CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/vector3D.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector3D.cpp > CMakeFiles/CMU462.dir/vector3D.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/vector3D.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector3D.cpp -o CMakeFiles/CMU462.dir/vector3D.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o: ../CMU462/src/vector4D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/vector4D.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector4D.cpp

CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/vector4D.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector4D.cpp > CMakeFiles/CMU462.dir/vector4D.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/vector4D.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/vector4D.cpp -o CMakeFiles/CMU462.dir/vector4D.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o: ../CMU462/src/matrix3x3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/matrix3x3.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/matrix3x3.cpp

CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/matrix3x3.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/matrix3x3.cpp > CMakeFiles/CMU462.dir/matrix3x3.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/matrix3x3.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/matrix3x3.cpp -o CMakeFiles/CMU462.dir/matrix3x3.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o: ../CMU462/src/matrix4x4.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/matrix4x4.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/matrix4x4.cpp

CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/matrix4x4.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/matrix4x4.cpp > CMakeFiles/CMU462.dir/matrix4x4.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/matrix4x4.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/matrix4x4.cpp -o CMakeFiles/CMU462.dir/matrix4x4.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o: ../CMU462/src/quaternion.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/quaternion.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/quaternion.cpp

CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/quaternion.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/quaternion.cpp > CMakeFiles/CMU462.dir/quaternion.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/quaternion.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/quaternion.cpp -o CMakeFiles/CMU462.dir/quaternion.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o: ../CMU462/src/complex.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/complex.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/complex.cpp

CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/complex.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/complex.cpp > CMakeFiles/CMU462.dir/complex.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/complex.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/complex.cpp -o CMakeFiles/CMU462.dir/complex.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o: ../CMU462/src/color.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/color.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/color.cpp

CMU462/src/CMakeFiles/CMU462.dir/color.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/color.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/color.cpp > CMakeFiles/CMU462.dir/color.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/color.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/color.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/color.cpp -o CMakeFiles/CMU462.dir/color.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o: ../CMU462/src/spectrum.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/spectrum.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/spectrum.cpp

CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/spectrum.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/spectrum.cpp > CMakeFiles/CMU462.dir/spectrum.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/spectrum.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/spectrum.cpp -o CMakeFiles/CMU462.dir/spectrum.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o: ../CMU462/src/osdtext.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/osdtext.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/osdtext.cpp

CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/osdtext.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/osdtext.cpp > CMakeFiles/CMU462.dir/osdtext.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/osdtext.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/osdtext.cpp -o CMakeFiles/CMU462.dir/osdtext.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o: ../CMU462/src/osdfont.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/CMU462.dir/osdfont.c.o   -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/osdfont.c

CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/CMU462.dir/osdfont.c.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/osdfont.c > CMakeFiles/CMU462.dir/osdfont.c.i

CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/CMU462.dir/osdfont.c.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/osdfont.c -o CMakeFiles/CMU462.dir/osdfont.c.s

CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.requires

CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.provides: CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.provides

CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o


CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o: ../CMU462/src/viewer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/viewer.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/viewer.cpp

CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/viewer.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/viewer.cpp > CMakeFiles/CMU462.dir/viewer.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/viewer.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/viewer.cpp -o CMakeFiles/CMU462.dir/viewer.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o: ../CMU462/src/base64.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/base64.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/base64.cpp

CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/base64.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/base64.cpp > CMakeFiles/CMU462.dir/base64.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/base64.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/base64.cpp -o CMakeFiles/CMU462.dir/base64.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o: ../CMU462/src/lodepng.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/lodepng.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/lodepng.cpp

CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/lodepng.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/lodepng.cpp > CMakeFiles/CMU462.dir/lodepng.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/lodepng.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/lodepng.cpp -o CMakeFiles/CMU462.dir/lodepng.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o


CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o: CMU462/src/CMakeFiles/CMU462.dir/flags.make
CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o: ../CMU462/src/tinyxml2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CMU462.dir/tinyxml2.cpp.o -c /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/tinyxml2.cpp

CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CMU462.dir/tinyxml2.cpp.i"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/tinyxml2.cpp > CMakeFiles/CMU462.dir/tinyxml2.cpp.i

CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CMU462.dir/tinyxml2.cpp.s"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/zach/Documents/15462/asst3_pathtracer/CMU462/src/tinyxml2.cpp -o CMakeFiles/CMU462.dir/tinyxml2.cpp.s

CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.requires:

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.requires

CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.provides: CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.requires
	$(MAKE) -f CMU462/src/CMakeFiles/CMU462.dir/build.make CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.provides.build
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.provides

CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.provides.build: CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o


# Object files for target CMU462
CMU462_OBJECTS = \
"CMakeFiles/CMU462.dir/vector2D.cpp.o" \
"CMakeFiles/CMU462.dir/vector3D.cpp.o" \
"CMakeFiles/CMU462.dir/vector4D.cpp.o" \
"CMakeFiles/CMU462.dir/matrix3x3.cpp.o" \
"CMakeFiles/CMU462.dir/matrix4x4.cpp.o" \
"CMakeFiles/CMU462.dir/quaternion.cpp.o" \
"CMakeFiles/CMU462.dir/complex.cpp.o" \
"CMakeFiles/CMU462.dir/color.cpp.o" \
"CMakeFiles/CMU462.dir/spectrum.cpp.o" \
"CMakeFiles/CMU462.dir/osdtext.cpp.o" \
"CMakeFiles/CMU462.dir/osdfont.c.o" \
"CMakeFiles/CMU462.dir/viewer.cpp.o" \
"CMakeFiles/CMU462.dir/base64.cpp.o" \
"CMakeFiles/CMU462.dir/lodepng.cpp.o" \
"CMakeFiles/CMU462.dir/tinyxml2.cpp.o"

# External object files for target CMU462
CMU462_EXTERNAL_OBJECTS =

CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/build.make
CMU462/src/libCMU462.a: CMU462/src/CMakeFiles/CMU462.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zach/Documents/15462/asst3_pathtracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking CXX static library libCMU462.a"
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && $(CMAKE_COMMAND) -P CMakeFiles/CMU462.dir/cmake_clean_target.cmake
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CMU462.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMU462/src/CMakeFiles/CMU462.dir/build: CMU462/src/libCMU462.a

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/build

CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/vector2D.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/vector3D.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/vector4D.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/matrix3x3.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/matrix4x4.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/quaternion.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/complex.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/color.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/spectrum.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/osdtext.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/osdfont.c.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/viewer.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/base64.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/lodepng.cpp.o.requires
CMU462/src/CMakeFiles/CMU462.dir/requires: CMU462/src/CMakeFiles/CMU462.dir/tinyxml2.cpp.o.requires

.PHONY : CMU462/src/CMakeFiles/CMU462.dir/requires

CMU462/src/CMakeFiles/CMU462.dir/clean:
	cd /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src && $(CMAKE_COMMAND) -P CMakeFiles/CMU462.dir/cmake_clean.cmake
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/clean

CMU462/src/CMakeFiles/CMU462.dir/depend:
	cd /home/zach/Documents/15462/asst3_pathtracer/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zach/Documents/15462/asst3_pathtracer /home/zach/Documents/15462/asst3_pathtracer/CMU462/src /home/zach/Documents/15462/asst3_pathtracer/build /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src /home/zach/Documents/15462/asst3_pathtracer/build/CMU462/src/CMakeFiles/CMU462.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMU462/src/CMakeFiles/CMU462.dir/depend

