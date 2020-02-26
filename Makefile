# CARGS = --std=c++14 -g


MY_ARGS = -Wall -Werror
ARGS_A = -std=gnu++11 -Wl,-rpath,/home/ben/Documents/ncbi/GCC800-DebugMT64/lib -L. -Wl,--enable-new-dtags -Wl,-export-dynamic -pthread -g

OBJS = invoker.o util/util.o blast/blast.o prospector/prospector.o prospector/dlinked.o

LIBS = -L/home/ben/Documents/ncbi/GCC800-DebugMT64/lib -lblastinput-static -lncbi_xloader_blastdb_rmt-static -lncbi_xloader_blastdb-static -lxblastformat-static -lalign_format-static -ltaxon1-static -lblastdb_format-static -lgene_info-static -lxformat-static -lxcleanup-static -lgbseq-static -lmlacli-static -lmla-static -lmedlars-static -lpubmed-static -lvalid-static -lxobjedit-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -ltaxon3-static -lxalnmgr-static -lblastxml-static -lblastxml2-static -lxcgi-static -lxhtml-static -lproteinkmer-static -lxblast-static -lxalgoblastdbindex-static -lcomposition_adjustment-static -lxalgodustmask-static -lxalgowinmask-static -lseqmasks_io-static -lseqdb-static -lblast_services-static -lxalnmgr-static -lxobjutil-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -lxnetblastcli-static -lxnetblast-static -lblastdb-static -lscoremat-static -ltables-static -llmdb-static -lncbi_xloader_genbank-static -lncbi_xreader_id1-static -lncbi_xreader_id2-static -lncbi_xreader_cache-static -ldbapi_driver-static -lncbi_xreader-static -lxconnect-static -lid1-static -lid2-static -lxobjmgr-static -lgenome_collection-static -lseqedit-static -lseqsplit-static -lsubmit-static -lseqset-static -lseq-static -lseqcode-static -lsequtil-static -lpub-static -lmedline-static -lbiblio-static -lgeneral-static -lxser-static -lxutil-static -lxncbi-static -lxcompress-static -lz-static -lbz2-static -lpthread -lnsl -ldl -ldl -lm -lpthread

CUDA_LIB = -L/usr/local/cuda/lib64 -lcudart


# libraries.o:
	# ld -r -o libraries.o --whole-archive $(LIBS)

run: invoker.out
	./invoker.out

invoker.out: clean $(OBJS)
	g++ -g $(MY_ARGS) $(ARGS_A) $(OBJS) $(LIBS) $(CUDA_LIB) -o invoker.out -ftime-report

invoker.o: invoker.cpp util/util.cpp prospector/prospector.cu blast/blast.cpp
	cd util && make
	cd prospector && make
	cd blast && make
	g++ -c -g invoker.cpp -o invoker.o


clean:
	cd util && make clean
	cd prospector && make clean
	cd blast && make clean
	rm -f invoker.o
	rm -f invoker.out

	# rm $(OBJS)