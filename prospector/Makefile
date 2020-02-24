NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g
# OBJ_DIR = ../obj

lib: dlinked.o prospector.o
	
dlinked.o: prospector.o
	$(NVCC) -dlink -o dlinked.o prospector.o
	
prospector.o:
	$(NVCC) $(CARGS) -dc -o prospector.o prospector.cu

clean:
	rm -fv *.o *.out
