
CARGS = --std=c++14 -g

invoker.out: device_linked
	g++ $(CARGS) -o invoker.out invoker.cpp ../prospector/device_linked.o ../prospector/util.o ../prospector/prospector.o -L../ -L/usr/local/cuda/lib64 -lcudart

device_linked:
	cd ../prospector && make device_linked.o

clean:
	rm -fv *.o *.so *.out *.a
	cd ../prospector && make clean
