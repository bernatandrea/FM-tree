# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /cygdrive/c/Users/anelh/.CLion2019.2/system/cygwin_cmake/bin/cmake.exe

# The command to remove a file.
RM = /cygdrive/c/Users/anelh/.CLion2019.2/system/cygwin_cmake/bin/cmake.exe -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cygdrive/c/Users/anelh/CLionProjects/FMTree

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/FMtree.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FMtree.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FMtree.dir/flags.make

CMakeFiles/FMtree.dir/main.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/FMtree.dir/main.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/main.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/main.cpp

CMakeFiles/FMtree.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/main.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/main.cpp > CMakeFiles/FMtree.dir/main.cpp.i

CMakeFiles/FMtree.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/main.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/main.cpp -o CMakeFiles/FMtree.dir/main.cpp.s

CMakeFiles/FMtree.dir/divsufsort.c.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/divsufsort.c.o: ../divsufsort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/FMtree.dir/divsufsort.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/FMtree.dir/divsufsort.c.o   -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/divsufsort.c

CMakeFiles/FMtree.dir/divsufsort.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/FMtree.dir/divsufsort.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/divsufsort.c > CMakeFiles/FMtree.dir/divsufsort.c.i

CMakeFiles/FMtree.dir/divsufsort.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/FMtree.dir/divsufsort.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/divsufsort.c -o CMakeFiles/FMtree.dir/divsufsort.c.s

CMakeFiles/FMtree.dir/fm-tree.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/fm-tree.cpp.o: ../fm-tree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/FMtree.dir/fm-tree.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/fm-tree.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/fm-tree.cpp

CMakeFiles/FMtree.dir/fm-tree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/fm-tree.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/fm-tree.cpp > CMakeFiles/FMtree.dir/fm-tree.cpp.i

CMakeFiles/FMtree.dir/fm-tree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/fm-tree.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/fm-tree.cpp -o CMakeFiles/FMtree.dir/fm-tree.cpp.s

CMakeFiles/FMtree.dir/test_me.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/test_me.cpp.o: ../test_me.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/FMtree.dir/test_me.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/test_me.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/test_me.cpp

CMakeFiles/FMtree.dir/test_me.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/test_me.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/test_me.cpp > CMakeFiles/FMtree.dir/test_me.cpp.i

CMakeFiles/FMtree.dir/test_me.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/test_me.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/test_me.cpp -o CMakeFiles/FMtree.dir/test_me.cpp.s

CMakeFiles/FMtree.dir/lib/bitset.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/bitset.cpp.o: ../lib/bitset.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/FMtree.dir/lib/bitset.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/bitset.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitset.cpp

CMakeFiles/FMtree.dir/lib/bitset.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/bitset.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitset.cpp > CMakeFiles/FMtree.dir/lib/bitset.cpp.i

CMakeFiles/FMtree.dir/lib/bitset.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/bitset.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitset.cpp -o CMakeFiles/FMtree.dir/lib/bitset.cpp.s

CMakeFiles/FMtree.dir/lib/bitmask.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/bitmask.cpp.o: ../lib/bitmask.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/FMtree.dir/lib/bitmask.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/bitmask.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask.cpp

CMakeFiles/FMtree.dir/lib/bitmask.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/bitmask.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask.cpp > CMakeFiles/FMtree.dir/lib/bitmask.cpp.i

CMakeFiles/FMtree.dir/lib/bitmask.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/bitmask.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask.cpp -o CMakeFiles/FMtree.dir/lib/bitmask.cpp.s

CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.o: ../lib/bitmask_vector.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask_vector.cpp

CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask_vector.cpp > CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.i

CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask_vector.cpp -o CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.s

CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.o: ../lib/bitmask_bitset.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask_bitset.cpp

CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask_bitset.cpp > CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.i

CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/bitmask_bitset.cpp -o CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.s

CMakeFiles/FMtree.dir/lib/wavelet.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/wavelet.cpp.o: ../lib/wavelet.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/FMtree.dir/lib/wavelet.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/wavelet.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/wavelet.cpp

CMakeFiles/FMtree.dir/lib/wavelet.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/wavelet.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/wavelet.cpp > CMakeFiles/FMtree.dir/lib/wavelet.cpp.i

CMakeFiles/FMtree.dir/lib/wavelet.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/wavelet.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/wavelet.cpp -o CMakeFiles/FMtree.dir/lib/wavelet.cpp.s

CMakeFiles/FMtree.dir/lib/rb_node.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/rb_node.cpp.o: ../lib/rb_node.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/FMtree.dir/lib/rb_node.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/rb_node.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/rb_node.cpp

CMakeFiles/FMtree.dir/lib/rb_node.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/rb_node.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/rb_node.cpp > CMakeFiles/FMtree.dir/lib/rb_node.cpp.i

CMakeFiles/FMtree.dir/lib/rb_node.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/rb_node.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/rb_node.cpp -o CMakeFiles/FMtree.dir/lib/rb_node.cpp.s

CMakeFiles/FMtree.dir/lib/rb_tree.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/rb_tree.cpp.o: ../lib/rb_tree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/FMtree.dir/lib/rb_tree.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/rb_tree.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/rb_tree.cpp

CMakeFiles/FMtree.dir/lib/rb_tree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/rb_tree.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/rb_tree.cpp > CMakeFiles/FMtree.dir/lib/rb_tree.cpp.i

CMakeFiles/FMtree.dir/lib/rb_tree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/rb_tree.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/rb_tree.cpp -o CMakeFiles/FMtree.dir/lib/rb_tree.cpp.s

CMakeFiles/FMtree.dir/lib/lookup_list.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/lookup_list.cpp.o: ../lib/lookup_list.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/FMtree.dir/lib/lookup_list.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/lookup_list.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/lookup_list.cpp

CMakeFiles/FMtree.dir/lib/lookup_list.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/lookup_list.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/lookup_list.cpp > CMakeFiles/FMtree.dir/lib/lookup_list.cpp.i

CMakeFiles/FMtree.dir/lib/lookup_list.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/lookup_list.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/lookup_list.cpp -o CMakeFiles/FMtree.dir/lib/lookup_list.cpp.s

CMakeFiles/FMtree.dir/lib/data.cpp.o: CMakeFiles/FMtree.dir/flags.make
CMakeFiles/FMtree.dir/lib/data.cpp.o: ../lib/data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/FMtree.dir/lib/data.cpp.o"
	/usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FMtree.dir/lib/data.cpp.o -c /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/data.cpp

CMakeFiles/FMtree.dir/lib/data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FMtree.dir/lib/data.cpp.i"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/data.cpp > CMakeFiles/FMtree.dir/lib/data.cpp.i

CMakeFiles/FMtree.dir/lib/data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FMtree.dir/lib/data.cpp.s"
	/usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/anelh/CLionProjects/FMTree/lib/data.cpp -o CMakeFiles/FMtree.dir/lib/data.cpp.s

# Object files for target FMtree
FMtree_OBJECTS = \
"CMakeFiles/FMtree.dir/main.cpp.o" \
"CMakeFiles/FMtree.dir/divsufsort.c.o" \
"CMakeFiles/FMtree.dir/fm-tree.cpp.o" \
"CMakeFiles/FMtree.dir/test_me.cpp.o" \
"CMakeFiles/FMtree.dir/lib/bitset.cpp.o" \
"CMakeFiles/FMtree.dir/lib/bitmask.cpp.o" \
"CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.o" \
"CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.o" \
"CMakeFiles/FMtree.dir/lib/wavelet.cpp.o" \
"CMakeFiles/FMtree.dir/lib/rb_node.cpp.o" \
"CMakeFiles/FMtree.dir/lib/rb_tree.cpp.o" \
"CMakeFiles/FMtree.dir/lib/lookup_list.cpp.o" \
"CMakeFiles/FMtree.dir/lib/data.cpp.o"

# External object files for target FMtree
FMtree_EXTERNAL_OBJECTS =

FMtree.exe: CMakeFiles/FMtree.dir/main.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/divsufsort.c.o
FMtree.exe: CMakeFiles/FMtree.dir/fm-tree.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/test_me.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/bitset.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/bitmask.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/bitmask_vector.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/bitmask_bitset.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/wavelet.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/rb_node.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/rb_tree.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/lookup_list.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/lib/data.cpp.o
FMtree.exe: CMakeFiles/FMtree.dir/build.make
FMtree.exe: CMakeFiles/FMtree.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX executable FMtree.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FMtree.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FMtree.dir/build: FMtree.exe

.PHONY : CMakeFiles/FMtree.dir/build

CMakeFiles/FMtree.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FMtree.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FMtree.dir/clean

CMakeFiles/FMtree.dir/depend:
	cd /cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cygdrive/c/Users/anelh/CLionProjects/FMTree /cygdrive/c/Users/anelh/CLionProjects/FMTree /cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug /cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug /cygdrive/c/Users/anelh/CLionProjects/FMTree/cmake-build-debug/CMakeFiles/FMtree.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FMtree.dir/depend

