CC = g++

# For program run

CC_OPTION = -fomit-frame-pointer -W -Wall -Winline -DNDBUG -O3 -mpopcnt -lpsapi

#CC_OPTION = -W -Wall -Winline -DNDBUG -O3 -mpopcnt -pg

# For debug
# CC_OPTION = -g -fomit-frame-pointer -W -Wall -Winline


ALL_OPS = tryall core Makeours Clear
ALL:$(ALL_OPS)

core:
	$(CC) -c lib/bitset.cpp lib/bitmask.cpp lib/bitmask_vector.cpp lib/bitmask_bitset.cpp lib/wavelet.cpp lib/rb_node.cpp lib/rb_tree.cpp lib/lookup_list.cpp lib/data.cpp divsufsort.c fm-tree.cpp $(CC_OPTION)
tryall:
	$(CC) -c main.cpp $(CC_OPTION)
Makeours:
	$(CC) -o FMtree main.o bitset.o bitmask.o bitmask_vector.o bitmask_bitset.o wavelet.o rb_node.o rb_tree.o lookup_list.o data.o divsufsort.o fm-tree.o $(CC_OPTION)
Clear:
	rm *.o
     
      
