# cmake_minimum_required(VERSION 3.16 FATAL_ERROR)



# project(crispr LANGUAGES CXX CUDA)

# message("Your C++ compiler supports these C++ features:")
# foreach(i ${CMAKE_CXX_COMPILE_FEATURES})
#   message("${i}")
# endforeach()


# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)


# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fuse-ld=lld")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fuse-ld=lld")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()



# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G -g -Xcompiler -fopenmp")
# set(CMAKE_BUILD_TYPE Debug)



# include_directories(
#         /usr/local/cuda/include
#         /home/ben/Documents/ncbi/include
#         /home/ben/Documents/ncbi/GCC800-DebugMT64/inc
#         /home/ben/Documents/ncbi/GCC800-DebugMT64/lib
# )


# add_library(
#         prospector STATIC
#         prospector/prospector.cu
#         prospector/prospector.h
# )


# add_library(
#         blast STATIC
#         blast/blast.cpp
# )


# target_include_directories(
#         blast
#         PRIVATE
#         /home/ben/Documents/ncbi/GCC800-DebugMT64/inc
#         /home/ben/Documents/ncbi/include
#         /home/ben/Documents/ncbi/include/internal
# )


# # target_compile_features(prospector PUBLIC cxx_std_11)
# set_target_properties(prospector PROPERTIES CUDA_SEPARABLE_COMPILATION ON)




# add_executable(invoker invoker.cpp)

# target_link_directories(
#         invoker
#         PRIVATE
#         /home/ben/Documents/ncbi/GCC800-DebugMT64/lib
# )


# set_property(TARGET invoker PROPERTY CUDA_SEPARABLE_COMPILATION ON)


# target_link_libraries(
#         invoker
#         PRIVATE
        
        
#         prospector
#         # util
#         blast

#         blastinput-static
#         ncbi_xloader_blastdb_rmt-static
#         ncbi_xloader_blastdb-static
#         xblastformat-static
#         align_format-static
#         taxon1-static
#         blastdb_format-static
#         gene_info-static
#         xformat-static
#         xcleanup-static
#         gbseq-static
#         mlacli-static
#         mla-static
#         medlars-static
#         pubmed-static
#         valid-static
#         xobjedit-static
#         xobjread-static
#         variation-static
#         creaders-static
#         submit-static
#         taxon3-static
#         xalnmgr-static
#         blastxml-static
#         blastxml2-static
#         xcgi-static
#         xhtml-static
#         proteinkmer-static
#         xblast-static
#         xalgoblastdbindex-static
#         composition_adjustment-static 
#         xalgodustmask-static
#         xalgowinmask-static
#         seqmasks_io-static 
#         seqdb-static
#         blast_services-static
#         xalnmgr-static
#         xobjutil-static
#         xobjread-static
#         variation-static
#         creaders-static
#         submit-static
#         xnetblastcli-static
#         xnetblast-static
#         blastdb-static
#         scoremat-static
#         tables-static
#         lmdb-static
#         ncbi_xloader_genbank-static
#         ncbi_xreader_id1-static
#         ncbi_xreader_id2-static
#         ncbi_xreader_cache-static
#         dbapi_driver-static
#         ncbi_xreader-static
#         xconnect-static
#         id1-static
#         id2-static
#         xobjmgr-static
#         genome_collection-static
#         seqedit-static
#         seqsplit-static
#         submit-static
#         seqset-static
#         seq-static
#         seqcode-static
#         sequtil-static
#         pub-static
#         medline-static
#         biblio-static
#         general-static
#         xser-static
#         xutil-static
#         xncbi-static
#         xcompress-static
#         z-static
#         bz2-static
#         pthread
#         nsl
#         dl
#         dl
#         m
#         pthread
# )

