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
CMAKE_COMMAND = /home/erick/.local/lib/python3.8/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/erick/.local/lib/python3.8/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/erick/Desktop/a24_ift3355_tp2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/erick/Desktop/a24_ift3355_tp2/build

# Include any dependencies generated for this target.
include CMakeFiles/RAY.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/RAY.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/RAY.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/RAY.dir/flags.make

CMakeFiles/RAY.dir/src/main.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/main.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/main.cpp
CMakeFiles/RAY.dir/src/main.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/RAY.dir/src/main.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/main.cpp.o -MF CMakeFiles/RAY.dir/src/main.cpp.o.d -o CMakeFiles/RAY.dir/src/main.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/main.cpp

CMakeFiles/RAY.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/main.cpp > CMakeFiles/RAY.dir/src/main.cpp.i

CMakeFiles/RAY.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/main.cpp -o CMakeFiles/RAY.dir/src/main.cpp.s

CMakeFiles/RAY.dir/src/object.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/object.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/object.cpp
CMakeFiles/RAY.dir/src/object.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/RAY.dir/src/object.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/object.cpp.o -MF CMakeFiles/RAY.dir/src/object.cpp.o.d -o CMakeFiles/RAY.dir/src/object.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/object.cpp

CMakeFiles/RAY.dir/src/object.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/object.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/object.cpp > CMakeFiles/RAY.dir/src/object.cpp.i

CMakeFiles/RAY.dir/src/object.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/object.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/object.cpp -o CMakeFiles/RAY.dir/src/object.cpp.s

CMakeFiles/RAY.dir/src/parser.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/parser.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/parser.cpp
CMakeFiles/RAY.dir/src/parser.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/RAY.dir/src/parser.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/parser.cpp.o -MF CMakeFiles/RAY.dir/src/parser.cpp.o.d -o CMakeFiles/RAY.dir/src/parser.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/parser.cpp

CMakeFiles/RAY.dir/src/parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/parser.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/parser.cpp > CMakeFiles/RAY.dir/src/parser.cpp.i

CMakeFiles/RAY.dir/src/parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/parser.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/parser.cpp -o CMakeFiles/RAY.dir/src/parser.cpp.s

CMakeFiles/RAY.dir/src/raytracer.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/raytracer.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/raytracer.cpp
CMakeFiles/RAY.dir/src/raytracer.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/RAY.dir/src/raytracer.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/raytracer.cpp.o -MF CMakeFiles/RAY.dir/src/raytracer.cpp.o.d -o CMakeFiles/RAY.dir/src/raytracer.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/raytracer.cpp

CMakeFiles/RAY.dir/src/raytracer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/raytracer.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/raytracer.cpp > CMakeFiles/RAY.dir/src/raytracer.cpp.i

CMakeFiles/RAY.dir/src/raytracer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/raytracer.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/raytracer.cpp -o CMakeFiles/RAY.dir/src/raytracer.cpp.s

CMakeFiles/RAY.dir/src/container.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/container.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/container.cpp
CMakeFiles/RAY.dir/src/container.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/RAY.dir/src/container.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/container.cpp.o -MF CMakeFiles/RAY.dir/src/container.cpp.o.d -o CMakeFiles/RAY.dir/src/container.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/container.cpp

CMakeFiles/RAY.dir/src/container.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/container.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/container.cpp > CMakeFiles/RAY.dir/src/container.cpp.i

CMakeFiles/RAY.dir/src/container.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/container.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/container.cpp -o CMakeFiles/RAY.dir/src/container.cpp.s

CMakeFiles/RAY.dir/src/aabb.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/aabb.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/aabb.cpp
CMakeFiles/RAY.dir/src/aabb.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/RAY.dir/src/aabb.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/aabb.cpp.o -MF CMakeFiles/RAY.dir/src/aabb.cpp.o.d -o CMakeFiles/RAY.dir/src/aabb.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/aabb.cpp

CMakeFiles/RAY.dir/src/aabb.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/aabb.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/aabb.cpp > CMakeFiles/RAY.dir/src/aabb.cpp.i

CMakeFiles/RAY.dir/src/aabb.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/aabb.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/aabb.cpp -o CMakeFiles/RAY.dir/src/aabb.cpp.s

CMakeFiles/RAY.dir/src/resource_manager.cpp.o: CMakeFiles/RAY.dir/flags.make
CMakeFiles/RAY.dir/src/resource_manager.cpp.o: /home/erick/Desktop/a24_ift3355_tp2/src/resource_manager.cpp
CMakeFiles/RAY.dir/src/resource_manager.cpp.o: CMakeFiles/RAY.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/RAY.dir/src/resource_manager.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/RAY.dir/src/resource_manager.cpp.o -MF CMakeFiles/RAY.dir/src/resource_manager.cpp.o.d -o CMakeFiles/RAY.dir/src/resource_manager.cpp.o -c /home/erick/Desktop/a24_ift3355_tp2/src/resource_manager.cpp

CMakeFiles/RAY.dir/src/resource_manager.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/RAY.dir/src/resource_manager.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erick/Desktop/a24_ift3355_tp2/src/resource_manager.cpp > CMakeFiles/RAY.dir/src/resource_manager.cpp.i

CMakeFiles/RAY.dir/src/resource_manager.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/RAY.dir/src/resource_manager.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erick/Desktop/a24_ift3355_tp2/src/resource_manager.cpp -o CMakeFiles/RAY.dir/src/resource_manager.cpp.s

# Object files for target RAY
RAY_OBJECTS = \
"CMakeFiles/RAY.dir/src/main.cpp.o" \
"CMakeFiles/RAY.dir/src/object.cpp.o" \
"CMakeFiles/RAY.dir/src/parser.cpp.o" \
"CMakeFiles/RAY.dir/src/raytracer.cpp.o" \
"CMakeFiles/RAY.dir/src/container.cpp.o" \
"CMakeFiles/RAY.dir/src/aabb.cpp.o" \
"CMakeFiles/RAY.dir/src/resource_manager.cpp.o"

# External object files for target RAY
RAY_EXTERNAL_OBJECTS =

RAY: CMakeFiles/RAY.dir/src/main.cpp.o
RAY: CMakeFiles/RAY.dir/src/object.cpp.o
RAY: CMakeFiles/RAY.dir/src/parser.cpp.o
RAY: CMakeFiles/RAY.dir/src/raytracer.cpp.o
RAY: CMakeFiles/RAY.dir/src/container.cpp.o
RAY: CMakeFiles/RAY.dir/src/aabb.cpp.o
RAY: CMakeFiles/RAY.dir/src/resource_manager.cpp.o
RAY: CMakeFiles/RAY.dir/build.make
RAY: CMakeFiles/RAY.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable RAY"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/RAY.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/RAY.dir/build: RAY
.PHONY : CMakeFiles/RAY.dir/build

CMakeFiles/RAY.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/RAY.dir/cmake_clean.cmake
.PHONY : CMakeFiles/RAY.dir/clean

CMakeFiles/RAY.dir/depend:
	cd /home/erick/Desktop/a24_ift3355_tp2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erick/Desktop/a24_ift3355_tp2 /home/erick/Desktop/a24_ift3355_tp2 /home/erick/Desktop/a24_ift3355_tp2/build /home/erick/Desktop/a24_ift3355_tp2/build /home/erick/Desktop/a24_ift3355_tp2/build/CMakeFiles/RAY.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/RAY.dir/depend

