CARGS = --std=c++14 -g



ARGS_A = -std=gnu++11 -Wl,-rpath,/home/ben/Documents/ncbi/GCC800-DebugMT64/lib -L. -Wl,--enable-new-dtags -Wl,-export-dynamic -pthread -g


OBJS = util/util.o blast_old/blast_demo.o prospector/prospector.o prospector/dlinked.o

LIBS = -L/home/ben/Documents/ncbi/GCC800-DebugMT64/lib -lblastinput-static -lncbi_xloader_blastdb_rmt-static -lncbi_xloader_blastdb-static -lxblastformat-static -lalign_format-static -ltaxon1-static -lblastdb_format-static -lgene_info-static -lxformat-static -lxcleanup-static -lgbseq-static -lmlacli-static -lmla-static -lmedlars-static -lpubmed-static -lvalid-static -lxobjedit-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -ltaxon3-static -lxalnmgr-static -lblastxml-static -lblastxml2-static -lxcgi-static -lxhtml-static -lproteinkmer-static -lxblast-static -lxalgoblastdbindex-static -lcomposition_adjustment-static -lxalgodustmask-static -lxalgowinmask-static -lseqmasks_io-static -lseqdb-static -lblast_services-static -lxalnmgr-static -lxobjutil-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -lxnetblastcli-static -lxnetblast-static -lblastdb-static -lscoremat-static -ltables-static -llmdb-static -lncbi_xloader_genbank-static -lncbi_xreader_id1-static -lncbi_xreader_id2-static -lncbi_xreader_cache-static -ldbapi_driver-static -lncbi_xreader-static -lxconnect-static -lid1-static -lid2-static -lxobjmgr-static -lgenome_collection-static -lseqedit-static -lseqsplit-static -lsubmit-static -lseqset-static -lseq-static -lseqcode-static -lsequtil-static -lpub-static -lmedline-static -lbiblio-static -lgeneral-static -lxser-static -lxutil-static -lxncbi-static -lxcompress-static -lz-static -lbz2-static -lpthread -lnsl -ldl -ldl -lm -lpthread

CUDA_LIB = -L/usr/local/cuda/lib64 -lcudart

# g++ $(CARGS) -o invoker.out invoker.cpp blast_old/blast_demo.o


invoker.o: blast_old/blast_demo.cpp prospector/prospector.cu util/util.cpp
	/usr/bin/g++ $(ARGS_A) $(OBJS) $(LIBS) $(CUDA_LIB) invoker.cpp -o invoker.out




blast.o: blast/blast_demo.cpp util/util.cpp
	cd blast && make


prospector.o: prospector/prospector.cu util/util.cpp
	cd prospector && make

