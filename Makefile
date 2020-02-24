CARGS = --std=c++14 -g



invoker.o: blast/blast_demo.cpp prospector/prospector.cu util/util.cpp
	g++ $(CARGS) -o invoker.cpp blast/blast_demo.o



blast.o: blast/blast_demo.cpp util/util.cpp
	cd blast && make


prospector.o: prospector/prospector.cu util/util.cpp
	cd prospector && make

