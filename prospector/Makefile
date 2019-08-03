NVCC = /usr/local/cuda-10.1/bin/nvcc


build:
	$(NVCC) --version
	$(NVCC) --std=c++14 -g prospector.cu Sequence.cpp  -I . -o prospector.o
