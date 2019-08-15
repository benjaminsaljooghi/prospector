NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g

.PHONY: version clean
binaries = program.out

program: clean version util.o prospector.o
	$(NVCC) $(CARGS) prospector.o util.o -o program.out

version:
	$(NVCC) --version

util.o:
	$(NVCC) $(CARGS) --device-c util.cu -o util.o

prospector.o:
	$(NVCC) $(CARGS) --device-c prospector.cu -o prospector.o

clean:
	rm -fv $(binaries) *.o