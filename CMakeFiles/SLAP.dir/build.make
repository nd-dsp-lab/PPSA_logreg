# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /home/jzhao7/AccuracyPSA/PSAapp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jzhao7/AccuracyPSA/PSAapp/build

# Include any dependencies generated for this target.
include CMakeFiles/SLAP.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/SLAP.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/SLAP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SLAP.dir/flags.make

CMakeFiles/SLAP.dir/main.cpp.o: CMakeFiles/SLAP.dir/flags.make
CMakeFiles/SLAP.dir/main.cpp.o: ../main.cpp
CMakeFiles/SLAP.dir/main.cpp.o: CMakeFiles/SLAP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SLAP.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SLAP.dir/main.cpp.o -MF CMakeFiles/SLAP.dir/main.cpp.o.d -o CMakeFiles/SLAP.dir/main.cpp.o -c /home/jzhao7/AccuracyPSA/PSAapp/main.cpp

CMakeFiles/SLAP.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SLAP.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jzhao7/AccuracyPSA/PSAapp/main.cpp > CMakeFiles/SLAP.dir/main.cpp.i

CMakeFiles/SLAP.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SLAP.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jzhao7/AccuracyPSA/PSAapp/main.cpp -o CMakeFiles/SLAP.dir/main.cpp.s

CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o: CMakeFiles/SLAP.dir/flags.make
CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o: ../slaprns-scheme.cpp
CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o: CMakeFiles/SLAP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o -MF CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o.d -o CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o -c /home/jzhao7/AccuracyPSA/PSAapp/slaprns-scheme.cpp

CMakeFiles/SLAP.dir/slaprns-scheme.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SLAP.dir/slaprns-scheme.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jzhao7/AccuracyPSA/PSAapp/slaprns-scheme.cpp > CMakeFiles/SLAP.dir/slaprns-scheme.cpp.i

CMakeFiles/SLAP.dir/slaprns-scheme.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SLAP.dir/slaprns-scheme.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jzhao7/AccuracyPSA/PSAapp/slaprns-scheme.cpp -o CMakeFiles/SLAP.dir/slaprns-scheme.cpp.s

CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o: CMakeFiles/SLAP.dir/flags.make
CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o: ../PSA-cryptocontext.cpp
CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o: CMakeFiles/SLAP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o -MF CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o.d -o CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o -c /home/jzhao7/AccuracyPSA/PSAapp/PSA-cryptocontext.cpp

CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jzhao7/AccuracyPSA/PSAapp/PSA-cryptocontext.cpp > CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.i

CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jzhao7/AccuracyPSA/PSAapp/PSA-cryptocontext.cpp -o CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.s

CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o: CMakeFiles/SLAP.dir/flags.make
CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o: ../PSA-base-scheme.cpp
CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o: CMakeFiles/SLAP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o -MF CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o.d -o CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o -c /home/jzhao7/AccuracyPSA/PSAapp/PSA-base-scheme.cpp

CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jzhao7/AccuracyPSA/PSAapp/PSA-base-scheme.cpp > CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.i

CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jzhao7/AccuracyPSA/PSAapp/PSA-base-scheme.cpp -o CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.s

CMakeFiles/SLAP.dir/utils.cpp.o: CMakeFiles/SLAP.dir/flags.make
CMakeFiles/SLAP.dir/utils.cpp.o: ../utils.cpp
CMakeFiles/SLAP.dir/utils.cpp.o: CMakeFiles/SLAP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/SLAP.dir/utils.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SLAP.dir/utils.cpp.o -MF CMakeFiles/SLAP.dir/utils.cpp.o.d -o CMakeFiles/SLAP.dir/utils.cpp.o -c /home/jzhao7/AccuracyPSA/PSAapp/utils.cpp

CMakeFiles/SLAP.dir/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SLAP.dir/utils.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jzhao7/AccuracyPSA/PSAapp/utils.cpp > CMakeFiles/SLAP.dir/utils.cpp.i

CMakeFiles/SLAP.dir/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SLAP.dir/utils.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jzhao7/AccuracyPSA/PSAapp/utils.cpp -o CMakeFiles/SLAP.dir/utils.cpp.s

# Object files for target SLAP
SLAP_OBJECTS = \
"CMakeFiles/SLAP.dir/main.cpp.o" \
"CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o" \
"CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o" \
"CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o" \
"CMakeFiles/SLAP.dir/utils.cpp.o"

# External object files for target SLAP
SLAP_EXTERNAL_OBJECTS =

SLAP: CMakeFiles/SLAP.dir/main.cpp.o
SLAP: CMakeFiles/SLAP.dir/slaprns-scheme.cpp.o
SLAP: CMakeFiles/SLAP.dir/PSA-cryptocontext.cpp.o
SLAP: CMakeFiles/SLAP.dir/PSA-base-scheme.cpp.o
SLAP: CMakeFiles/SLAP.dir/utils.cpp.o
SLAP: CMakeFiles/SLAP.dir/build.make
SLAP: /usr/local/lib/libOPENFHEpke.so.1.2.0
SLAP: /usr/local/lib/libOPENFHEbinfhe.so.1.2.0
SLAP: /usr/local/lib/libOPENFHEcore.so.1.2.0
SLAP: CMakeFiles/SLAP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable SLAP"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SLAP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SLAP.dir/build: SLAP
.PHONY : CMakeFiles/SLAP.dir/build

CMakeFiles/SLAP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SLAP.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SLAP.dir/clean

CMakeFiles/SLAP.dir/depend:
	cd /home/jzhao7/AccuracyPSA/PSAapp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jzhao7/AccuracyPSA/PSAapp /home/jzhao7/AccuracyPSA/PSAapp /home/jzhao7/AccuracyPSA/PSAapp/build /home/jzhao7/AccuracyPSA/PSAapp/build /home/jzhao7/AccuracyPSA/PSAapp/build/CMakeFiles/SLAP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SLAP.dir/depend

