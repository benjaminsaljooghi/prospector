NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g
	
dlinked.o: prospector.o
	$(NVCC) -dlink -o dlinked.o prospector.o
	
prospector.o: prospector.cu
	$(NVCC) $(CARGS) -dc -o prospector.o prospector.cu

