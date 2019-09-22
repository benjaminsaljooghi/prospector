NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g

.PHONY: version clean



# program.out: device_linked.o
	# g++ -o program.out prospector.o util.o device_linked.o -L. -L/usr/local/cuda/lib64 -lcudart

program.out: dlinked.o prospector.o util.o
	g++ -o program.out prospector.o util.o dlinked.o -L. -L/usr/local/cuda/lib64 -lcudart

dlinked.o: prospector.o
	$(NVCC) -dlink -o dlinked.o prospector.o
	
prospector.o:
	$(NVCC) $(CARGS) -dc -o prospector.o prospector.cu

util.o:
	$(NVCC) $(CARGS) -c -o util.o util.cpp

version:
	$(NVCC) --version

clean:
	rm -fv *.o *.so *.out *.a