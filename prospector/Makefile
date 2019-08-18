NVCC = /usr/local/cuda/bin/nvcc
CARGS = --std=c++14 -g

.PHONY: version clean

# main.out: libgpu.a
	# g++ -o main.out main.cpp -L. -lgpu -L/usr/local/cuda/lib64 -lcudart	

device_linked.o: prospector.o util.o
	$(NVCC) -dlink -o device_linked.o prospector.o util.o
	
prospector.o:
	$(NVCC) $(CARGS) -dc -o prospector.o prospector.cu

util.o:
	$(NVCC) $(CARGS) -dc -o util.o util.cu

version:
	$(NVCC) --version

clean:
	rm -fv *.o *.so *.out *.a