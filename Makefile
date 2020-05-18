OPT = -O3
CPP = g++-9 $(OPT) --std=gnu++2a -g -Wall -fopenmp
NVCC = /usr/local/cuda/bin/nvcc --std=c++14 -g -G -Xcompiler -fopenmp $(OPT) 
B = build

.PHONY: run build rebuild clean
run: build
	$(B)/main.out
build: $(B)/main.out
rebuild: clean build
clean:
	rm -fv $(B)/*.o $(B)/*.out

$(B)/util.o: util.*
	$(CPP) -c util.cpp -o $(B)/util.o

$(B)/crispr.o: crispr.* $(B)/util.o
	$(CPP) -c crispr.cpp -o $(B)/crispr.o

$(B)/cas.o: cas.* $(B)/util.o $(B)/crispr.o
	$(CPP) -c cas.cpp -o $(B)/cas.o

$(B)/dlinked.o: $(B)/prospector.o
	$(NVCC) -dlink $(B)/prospector.o -o $(B)/dlinked.o
$(B)/prospector.o: prospector.*
	$(NVCC) -dc prospector.cu -o $(B)/prospector.o



PROSP = $(B)/prospector.o $(B)/dlinked.o
LIB = -L/usr/local/cuda/lib64 -lcudart -L/usr/local/lib -lfmt

$(B)/main.out: main.* $(PROSP) $(B)/crispr.o $(B)/util.o $(B)/cas.o
	$(CPP) -c main.cpp -o $(B)/main.o
	$(CPP) $(B)/*.o $(LIB) -o $(B)/main.out -fuse-ld=lld

