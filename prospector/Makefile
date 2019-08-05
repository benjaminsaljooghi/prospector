NVCC = /usr/local/cuda-10.1/bin/nvcc

.PHONY: version clean
binaries = program

program: version util.o prospector.o
	$(NVCC) --std=c++14 -g prospector.o util.o -o program

version:
	$(NVCC) --version

util.o:
	$(NVCC) --std=c++14 -g --device-c util.cu -o util.o

prospector.o:
	$(NVCC) --std=c++14 -g --device-c prospector.cu -o prospector.o

clean:
	rm -f $(binaries) *.o