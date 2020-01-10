cmake_minimum_required(VERSION 3.15)
project(FMtree)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES
        # bitset
        lib/bitset.cpp lib/bitset.h

        # bitmask
        lib/bitmask.cpp lib/bitmask.h
        lib/bitmask_bitset.cpp lib/bitmask_bitset.h

        # wavelet
        lib/wavelet.cpp lib/wavelet.h

        # lookup tree
        lib/rank_select.h

	 # rb tree
        lib/rb_node.cpp lib/rb_node.h
        lib/rb_tree.cpp lib/rb_tree.h

        # data
        lib/data.cpp lib/data.h
        )

add_executable(FMtree main.cpp divsufsort.c divsufsort.h fm-tree.cpp fm-tree.h ${SOURCE_FILES})
