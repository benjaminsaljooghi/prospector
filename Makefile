OPTIMIZATION = -O0

CPP = clang++
CARGS = --std=c++17 -g -Wall -fopenmp -g

NVCC = /usr/local/cuda/bin/nvcc
NVCCARGS = --std=c++14 -g -G -Xcompiler -fopenmp

NCBI = /home/ben/lib/ncbi
INC_NCBI = -I. -I$(NCBI)/GCC800-DebugMT64/inc -I$(NCBI)/include -I$(NCBI)/include/internal

LIB_ARGS = -Wl,-rpath,$(NCBI)/GCC800-DebugMT64/lib -L. -Wl,--enable-new-dtags -Wl,-export-dynamic -pthread -g
LIB_NCBI = -L$(NCBI)/GCC800-DebugMT64/lib -lblastinput-static -lncbi_xloader_blastdb_rmt-static -lncbi_xloader_blastdb-static -lxblastformat-static -lalign_format-static -ltaxon1-static -lblastdb_format-static -lgene_info-static -lxformat-static -lxcleanup-static -lgbseq-static -lmlacli-static -lmla-static -lmedlars-static -lpubmed-static -lvalid-static -lxobjedit-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -ltaxon3-static -lxalnmgr-static -lblastxml-static -lblastxml2-static -lxcgi-static -lxhtml-static -lproteinkmer-static -lxblast-static -lxalgoblastdbindex-static -lcomposition_adjustment-static -lxalgodustmask-static -lxalgowinmask-static -lseqmasks_io-static -lseqdb-static -lblast_services-static -lxalnmgr-static -lxobjutil-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -lxnetblastcli-static -lxnetblast-static -lblastdb-static -lscoremat-static -ltables-static -llmdb-static -lncbi_xloader_genbank-static -lncbi_xreader_id1-static -lncbi_xreader_id2-static -lncbi_xreader_cache-static -ldbapi_driver-static -lncbi_xreader-static -lxconnect-static -lid1-static -lid2-static -lxobjmgr-static -lgenome_collection-static -lseqedit-static -lseqsplit-static -lsubmit-static -lseqset-static -lseq-static -lseqcode-static -lsequtil-static -lpub-static -lmedline-static -lbiblio-static -lgeneral-static -lxser-static -lxutil-static -lxncbi-static -lxcompress-static -lz-static -lbz2-static -lpthread -lnsl -ldl -ldl -lm -lpthread
LIB_CUDA = -L/usr/local/cuda/lib64 -lcudart

BLAST_ARGS = -c -Wno-format-y2k  -pthread -fPIC -D_DEBUG -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE   -D_MT -D_REENTRANT -D_THREAD_SAFE

BUILD = build
OBJS = $(BUILD)/invoker.o $(BUILD)/dlinked.o $(BUILD)/prospector.o $(BUILD)/blast.o
SRC = *.cpp


run: $(BUILD)/invoker.out
	./$(BUILD)/invoker.out


.PHONY: all_objs
all_objs: *.cpp *.h
	
	$(CPP) $(OPTIMIZATION) $(CARGS) -c util.cpp -o $(BUILD)/util.o
	$(CPP) $(OPTIMIZATION) $(CARGS) -c crispr.cpp -o $(BUILD)/crispr.o
	$(CPP) $(OPTIMZATION) $(CARGS) $(BLAST_ARGS) $(INC_NCBI) blast.cpp -o $(BUILD)/blast.o
	$(NVCC) $(NVCCARGS) -dc prospector.cu -o $(BUILD)/prospector.o
	$(NVCC) $(NVCCARGS) -dlink $(BUILD)/prospector.o -o $(BUILD)/dlinked.o
	$(CPP) $(OPTIMIZATION) $(CARGS) -c invoker.cpp -o $(BUILD)/invoker.o
	$(CPP) $(OPTIMIZATION) $(CARGS) $(BUILD)/*.o $(LIB_ARGS) $(LIB_NCBI) $(LIB_CUDA) -o invoker.out








.PHONY: clean
clean:
	rm -fv $(BUILD)/*.o $(BUILD)/*.out