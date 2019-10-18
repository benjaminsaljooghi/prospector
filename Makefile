


build: clean
	cd prospector && make
	cd blast && make



clean:
	cd prospector && make clean
	cd blast && make clean