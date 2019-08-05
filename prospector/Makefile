NVCC = /usr/local/cuda-10.1/bin/nvcc
ARGS = --std=c++14 -g

build:
	$(NVCC) --version
	$(NVCC) $(ARGS) -c util.cu -o util.o
	$(NVCC) $(ARGS) -c prospector.cu -o prospector.o
	$(NVCC) $(ARGS) prospector.o util.o -o program

	