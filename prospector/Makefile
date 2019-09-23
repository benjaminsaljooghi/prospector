NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g
OBJ_DIR = ../obj

.PHONY: version clean

program: dlinked prospector util
	g++ -o $(OBJ_DIR)/program.out $(OBJ_DIR)/prospector.o $(OBJ_DIR)/util.o $(OBJ_DIR)/dlinked.o -L. -L/usr/local/cuda/lib64 -lcudart

dlinked: prospector
	$(NVCC) -dlink -o $(OBJ_DIR)/dlinked.o $(OBJ_DIR)/prospector.o
	
prospector:
	$(NVCC) $(CARGS) -dc -o $(OBJ_DIR)/prospector.o prospector.cu

util:
	g++ $(CARGS) -c -o $(OBJ_DIR)/util.o util.cpp

clean:
	rm -fv $(OBJ_DIR)/*.o $(OBJ_DIR)/*.out