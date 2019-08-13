NVCC = /usr/local/cuda-10.1/bin/nvcc
CARGS = --std=c++14 -g

.PHONY: version clean
binaries = program

program: clean version util.o prospector.o
	$(NVCC) $(CARGS) prospector.o util.o -o program

version:
	$(NVCC) --version

util.o:
	$(NVCC) $(CARGS) --device-c util.cu -o util.o

prospector.o:
	$(NVCC) $(CARGS) --device-c prospector.cu -o prospector.o

clean:
	rm -f $(binaries) *.o