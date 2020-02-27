

MY_ARGS = -O3
ARGS = -std=gnu++11 -c  -Wall -Wno-format-y2k  -pthread -g -fPIC -D_DEBUG -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE   -D_MT -D_REENTRANT -D_THREAD_SAFE  
INC = -I. -I/home/ben/Documents/ncbi/GCC800-DebugMT64/inc -I/home/ben/Documents/ncbi/include -I/home/ben/Documents/ncbi/include/internal


blast.o: blast.cpp
	g++ $(MY_ARGS) $(ARGS) $(INC) blast.cpp -o blast.o 
