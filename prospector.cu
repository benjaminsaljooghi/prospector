#include "prospector.h"


#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"

#include <cassert>
#include <chrono>
#include <stdio.h>
#include "time.h"


#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define KERNEL_ARGS2(grid, block) <<< grid, block >>>
#define KERNEL_ARGS3(grid, block, sh_mem) <<< grid, block, sh_mem >>>
#define KERNEL_ARGS4(grid, block, sh_mem, stream) <<< grid, block, sh_mem, stream >>>
#else
#define CUDA_CALLABLE_MEMBER
#define KERNEL_ARGS2(grid, block)
#define KERNEL_ARGS3(grid, block, sh_mem)
#define KERNEL_ARGS4(grid, block, sh_mem, stream)
#endif

#define DEBUG 0

#define C_GRID 128
#define C_BLOCK 1024


cudaError_t checkCuda(cudaError_t result)
{
#if DEBUG == 1
    if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
    return result;
}

__device__ int __popc (ui x);

__device__ uc difference_gpu(const ui& _a, const ui& _b)
{
    ui _xor = (_a ^ _b);
    ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ui oddBits = _xor & 0x5555555555555555ull;
    ui comp = (evenBits >> 1) | oddBits;
    return __popc(comp);
}

__global__ void compute_qmap(const ui* genome_encoding, const ui genome_encoding_size, uc* qmap)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    for (ui query = thread_id; query < genome_encoding_size - 200; query += stride)
    {
        ui q = genome_encoding[query];
        for (ui i = 0; i < MAP_SIZE; i++)
        {
            ui t = genome_encoding[query + K_START + SPACER_SKIP + i];
            qmap[(query*MAP_SIZE) + i] = difference_gpu(q, t);
        }
    }
}


__device__ ui scheme(const char c)
{
    switch (c)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
    }
}

__global__ void compute_encoding(const char* genome, ui* genome_encoding, ui genome_size, ui genome_encoding_size)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    for (ui i = thread_id; i < genome_encoding_size; i += stride)
    {
        ui e = 0;
        for (int j = 0; j < SIZE; j++)
        {
            e |= scheme(genome[i + j]) << (j * BITS);
        }
        genome_encoding[i] = e;
    }
}



void Prospector::device_init()
{
    std::chrono::_V2::system_clock::time_point start = time();
    cudaFree(0);
    time(start, "device init");
}


Prospector::Encoding Prospector::get_genome_encoding(const char* genome, ui genome_size)
{
    std::chrono::_V2::system_clock::time_point start = time();

    cudaError er;
    
    char* d_genome;
    ui bytes_genome = sizeof(char) * genome_size;
    er = cudaMalloc(&d_genome, bytes_genome); checkCuda(er);
    er = cudaMemcpy(d_genome, genome, bytes_genome, cudaMemcpyHostToDevice); checkCuda(er);

    ui* d_encoding;
    ui encoding_size = genome_size - SIZE + 1;
    ui bytes_encoding = sizeof(ui) * encoding_size; 
    er = cudaMalloc(&d_encoding, bytes_encoding); checkCuda(er);

    compute_encoding KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_genome, d_encoding, genome_size, encoding_size);

    ui* encoding;
    er = cudaMallocHost(&encoding, bytes_encoding); checkCuda(er);

    cudaDeviceSynchronize();

    er = cudaMemcpy(encoding, d_encoding, bytes_encoding, cudaMemcpyDeviceToHost); checkCuda(er);

    time(start, "genome encoding total");

    return Encoding {
        encoding,
        d_encoding,
        encoding_size
    };
}


uc* Prospector::get_qmap(ui* d_encoding, ui encoding_size)
{
    assert(K_START >= SIZE);
    cudaError er; 
    auto start = time();

    ui bytes_qmap = sizeof(uc) * encoding_size * MAP_SIZE;
    
    uc* d_qmap;
    
    er = cudaMalloc(&d_qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap malloc");

    uc* qmap;
    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap mallochost");

    compute_qmap KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_encoding, encoding_size, d_qmap);

    cudaDeviceSynchronize();
    start = time(start, "kernel");

    er = cudaMemcpy(qmap, d_qmap, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er);
    start = time(start, "qmap memcpy");

    cudaFree(d_qmap);

    return qmap;
}



__global__ void compute_qmap3000(const ui* encoding, const ui encoding_size, const ui* queries, const ui queries_size, uc* qmap, const ui map_size)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    // for (ui query = thread_id; query < genome_encoding_size - 200; query += stride)
    for (ui query_i = thread_id; query_i < queries_size; query_i += stride)
    {
        const ui query = queries[query_i];
        const ui q = encoding[query];
        for (ui i = 0; i < map_size; i++)
        {
            const ui t = encoding[query + K_START + SPACER_SKIP + i];
            qmap[(query_i*map_size) + i] = difference_gpu(q, t);
        }
    }
}


uc* Prospector::get_qmap3000(const ui* d_encoding, const ui encoding_size, const ui* queries, const ui queries_size)
{
    ui map_size = 3000;

    ui bytes_queries = sizeof(uc) * queries_size;

    ui* d_queries;

    cudaMalloc(&d_queries, bytes_queries);
    cudaMemcpy(d_queries, queries, bytes_queries, cudaMemcpyHostToDevice);

    ui bytes_qmap = sizeof(uc) * queries_size * map_size;

    uc* d_qmap;
    uc* qmap;

    cudaMalloc(&d_qmap, bytes_qmap);
    cudaMallocHost(&qmap, bytes_qmap);

    compute_qmap3000 KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_encoding, encoding_size, queries, queries_size, d_qmap, map_size);
    cudaDeviceSynchronize();

    cudaMemcpy(qmap, d_qmap, bytes_qmap, cudaMemcpyDeviceToHost);
    cudaFree(d_qmap);

    return qmap;
}


