NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g

.PHONY: version clean

prospector.out: clean version util.o prospector.o
	$(NVCC) $(CARGS) prospector.o util.o -o prospector.out

prospector.o:
	$(NVCC) $(CARGS) --device-c prospector.cu -o prospector.o

util.o:
	$(NVCC) $(CARGS) --device-c util.cu -o util.o

version:
	$(NVCC) --version

clean:
	rm -fv *.o *.out