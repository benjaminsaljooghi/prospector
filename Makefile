
OPT = -O3

CPP = clang++ $(OPT) --std=gnu++2a -g -Wall -fopenmp
NVCC = /usr/local/cuda/bin/nvcc --std=c++14 -g -G -Xcompiler -fopenmp $(OPT) 

# NCBI = /home/ben/bin/ncbi
# INC_NCBI = -I. -I$(NCBI)/GCC800-DebugMT64/inc -I$(NCBI)/include -I$(NCBI)/include/internal

# LIB_ARGS = -Wl,-rpath,$(NCBI)/GCC800-DebugMT64/lib -L. -Wl,--enable-new-dtags -Wl,-export-dynamic -pthread -g
# LIB_NCBI = -L$(NCBI)/GCC800-DebugMT64/lib -lblastinput-static -lncbi_xloader_blastdb_rmt-static -lncbi_xloader_blastdb-static -lxblastformat-static -lalign_format-static -ltaxon1-static -lblastdb_format-static -lgene_info-static -lxformat-static -lxcleanup-static -lgbseq-static -lmlacli-static -lmla-static -lmedlars-static -lpubmed-static -lvalid-static -lxobjedit-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -ltaxon3-static -lxalnmgr-static -lblastxml-static -lblastxml2-static -lxcgi-static -lxhtml-static -lproteinkmer-static -lxblast-static -lxalgoblastdbindex-static -lcomposition_adjustment-static -lxalgodustmask-static -lxalgowinmask-static -lseqmasks_io-static -lseqdb-static -lblast_services-static -lxalnmgr-static -lxobjutil-static -lxobjread-static -lvariation-static -lcreaders-static -lsubmit-static -lxnetblastcli-static -lxnetblast-static -lblastdb-static -lscoremat-static -ltables-static -llmdb-static -lncbi_xloader_genbank-static -lncbi_xreader_id1-static -lncbi_xreader_id2-static -lncbi_xreader_cache-static -ldbapi_driver-static -lncbi_xreader-static -lxconnect-static -lid1-static -lid2-static -lxobjmgr-static -lgenome_collection-static -lseqedit-static -lseqsplit-static -lsubmit-static -lseqset-static -lseq-static -lseqcode-static -lsequtil-static -lpub-static -lmedline-static -lbiblio-static -lgeneral-static -lxser-static -lxutil-static -lxncbi-static -lxcompress-static -lz-static -lbz2-static -lpthread -lnsl -ldl -ldl -lm -lpthread
LIB_CUDA = -L/usr/local/cuda/lib64 -lcudart

# LIB = $(LIB_ARGS) $(LIB_NCBI) $(LIB_CUDA) -L/usr/local/lib -lfmt
# LIB = $(LIB_ARGS) $(LIB_NCBI) -L/usr/local/lib -lfmt
LIB_CPP = -L/usr/local/lib -lfmt

# BLAST_ARGS = -Wno-format-y2k  -pthread -fPIC -D_DEBUG -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE   -D_MT -D_REENTRANT -D_THREAD_SAFE

B = build

.PHONY: build rebuild clean

build: $(B)/main.out

rebuild: clean build

clean:
	rm -fv $(B)/*.o $(B)/*.out

$(B)/util.o: util.*
	$(CPP) -c util.cpp -o $(B)/util.o

$(B)/crispr.o: crispr.* $(B)/util.o
	$(CPP) -c crispr.cpp -o $(B)/crispr.o

# $(B)/blast.o: blast.* $(B)/util.o $(B)/crispr.o
	# $(CPP) -c $(BLAST_ARGS) $(INC_NCBI) blast.cpp -o $(B)/blast.o


# ALT GPU
PROSP = $(B)/prospector.o $(B)/dlinked.o
LIB_ALL = $(LIB_CPP) $(LIB_CUDA)
$(B)/dlinked.o: $(B)/prospector.o
	$(NVCC) -dlink $(B)/prospector.o -o $(B)/dlinked.o
$(B)/prospector.o: prospector.* $(B)/util.o $(B)/crispr.o
	$(NVCC) -dc prospector.cu -o $(B)/prospector.o

# ALT CPU
# PROSP = $(B)/prospector.o
# LIB_ALL = $(LIB_CPP)
# $(B)/prospector.o: prospector.* $(B)/util.o $(B)/crispr.o
# 	$(CPP) -c prospector.cpp -o $(B)/prospector.o

$(B)/main.out: main.* $(PROSP) $(B)/crispr.o $(B)/util.o
	$(CPP) -c main.cpp -o $(B)/main.o
	$(CPP) $(B)/*.o $(LIB_ALL) -o $(B)/main.out -fuse-ld=lld

