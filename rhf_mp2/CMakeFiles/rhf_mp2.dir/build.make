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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/khan/psi4-workspace/plugins/rhf_mp2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/khan/psi4-workspace/plugins/rhf_mp2

# Include any dependencies generated for this target.
include CMakeFiles/rhf_mp2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/rhf_mp2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rhf_mp2.dir/flags.make

CMakeFiles/rhf_mp2.dir/tensors.cc.o: CMakeFiles/rhf_mp2.dir/flags.make
CMakeFiles/rhf_mp2.dir/tensors.cc.o: tensors.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/khan/psi4-workspace/plugins/rhf_mp2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rhf_mp2.dir/tensors.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rhf_mp2.dir/tensors.cc.o -c /home/khan/psi4-workspace/plugins/rhf_mp2/tensors.cc

CMakeFiles/rhf_mp2.dir/tensors.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rhf_mp2.dir/tensors.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/khan/psi4-workspace/plugins/rhf_mp2/tensors.cc > CMakeFiles/rhf_mp2.dir/tensors.cc.i

CMakeFiles/rhf_mp2.dir/tensors.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rhf_mp2.dir/tensors.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/khan/psi4-workspace/plugins/rhf_mp2/tensors.cc -o CMakeFiles/rhf_mp2.dir/tensors.cc.s

CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.o: CMakeFiles/rhf_mp2.dir/flags.make
CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.o: rhf_mp2.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/khan/psi4-workspace/plugins/rhf_mp2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.o -c /home/khan/psi4-workspace/plugins/rhf_mp2/rhf_mp2.cc

CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/khan/psi4-workspace/plugins/rhf_mp2/rhf_mp2.cc > CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.i

CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/khan/psi4-workspace/plugins/rhf_mp2/rhf_mp2.cc -o CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.s

CMakeFiles/rhf_mp2.dir/main.cc.o: CMakeFiles/rhf_mp2.dir/flags.make
CMakeFiles/rhf_mp2.dir/main.cc.o: main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/khan/psi4-workspace/plugins/rhf_mp2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/rhf_mp2.dir/main.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rhf_mp2.dir/main.cc.o -c /home/khan/psi4-workspace/plugins/rhf_mp2/main.cc

CMakeFiles/rhf_mp2.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rhf_mp2.dir/main.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/khan/psi4-workspace/plugins/rhf_mp2/main.cc > CMakeFiles/rhf_mp2.dir/main.cc.i

CMakeFiles/rhf_mp2.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rhf_mp2.dir/main.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/khan/psi4-workspace/plugins/rhf_mp2/main.cc -o CMakeFiles/rhf_mp2.dir/main.cc.s

# Object files for target rhf_mp2
rhf_mp2_OBJECTS = \
"CMakeFiles/rhf_mp2.dir/tensors.cc.o" \
"CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.o" \
"CMakeFiles/rhf_mp2.dir/main.cc.o"

# External object files for target rhf_mp2
rhf_mp2_EXTERNAL_OBJECTS =

rhf_mp2.so: CMakeFiles/rhf_mp2.dir/tensors.cc.o
rhf_mp2.so: CMakeFiles/rhf_mp2.dir/rhf_mp2.cc.o
rhf_mp2.so: CMakeFiles/rhf_mp2.dir/main.cc.o
rhf_mp2.so: CMakeFiles/rhf_mp2.dir/build.make
rhf_mp2.so: /opt/psi4-1.3.2/lib/psi4/core.cpython-39-x86_64-linux-gnu.so
rhf_mp2.so: /usr/lib/gcc/x86_64-redhat-linux/10/libgomp.so
rhf_mp2.so: /usr/lib64/libpthread.so
rhf_mp2.so: CMakeFiles/rhf_mp2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/khan/psi4-workspace/plugins/rhf_mp2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared module rhf_mp2.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rhf_mp2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/rhf_mp2.dir/build: rhf_mp2.so

.PHONY : CMakeFiles/rhf_mp2.dir/build

CMakeFiles/rhf_mp2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rhf_mp2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rhf_mp2.dir/clean

CMakeFiles/rhf_mp2.dir/depend:
	cd /home/khan/psi4-workspace/plugins/rhf_mp2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/khan/psi4-workspace/plugins/rhf_mp2 /home/khan/psi4-workspace/plugins/rhf_mp2 /home/khan/psi4-workspace/plugins/rhf_mp2 /home/khan/psi4-workspace/plugins/rhf_mp2 /home/khan/psi4-workspace/plugins/rhf_mp2/CMakeFiles/rhf_mp2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rhf_mp2.dir/depend
