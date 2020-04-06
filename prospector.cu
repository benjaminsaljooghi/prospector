#include "prospector.h"


#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"

#include <cassert>
#include <chrono>
#include <stdio.h>


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

__device__ int __popc (unsigned int x);

__device__ unsigned char difference_gpu(const ui& _a, const ui& _b)
{
    ui _xor = (_a ^ _b);
    ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ui oddBits = _xor & 0x5555555555555555ull;
    ui comp = (evenBits >> 1) | oddBits;
    return __popc(comp);
}

__global__ void compute_qmap(const ui* genome_encoding, const ui genome_encoding_size, unsigned char* qmap)
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


std::chrono::_V2::system_clock::time_point time()
{
    return std::chrono::high_resolution_clock::now();
}

std::chrono::_V2::system_clock::time_point time(std::chrono::_V2::system_clock::time_point start, const char* message)
{
    auto curr = time();
    printf("%ldms %s\n", std::chrono::duration_cast<std::chrono::milliseconds>(curr - start).count(), message);
    return curr;
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

    ui* d_genome_encoding;
    ui genome_encoding_size = genome_size - SIZE + 1;
    ui bytes_genome_encoding = sizeof(ui) * genome_encoding_size; 
    er = cudaMalloc(&d_genome_encoding, bytes_genome_encoding); checkCuda(er);

    compute_encoding KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_genome, d_genome_encoding, genome_size, genome_encoding_size);

    ui* genome_encoding;
    er = cudaMallocHost(&genome_encoding, bytes_genome_encoding); checkCuda(er);

    cudaDeviceSynchronize();

    er = cudaMemcpy(genome_encoding, d_genome_encoding, bytes_genome_encoding, cudaMemcpyDeviceToHost); checkCuda(er);

    time(start, "genome encoding total");


    Prospector::Encoding encoding;

    encoding.encoding = genome_encoding;
    encoding.d_encoding = d_genome_encoding;

    return encoding;
}


uc* Prospector::get_qmap(ui* d_encoding, ui genome_encoding_size)
{
    assert(K_START >= SIZE);
    cudaError er; 
    std::chrono::_V2::system_clock::time_point start = time();


    ui bytes_qmap = sizeof(uc) * genome_encoding_size * MAP_SIZE;
    
    uc* d_qmap;
    
    er = cudaMalloc(&d_qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap malloc");

    uc* qmap;
    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap mallochost");

    compute_qmap KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_encoding, genome_encoding_size, d_qmap);

    cudaDeviceSynchronize();
    start = time(start, "kernel");

    er = cudaMemcpy(qmap, d_qmap, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er);
    start = time(start, "qmap memcpy");

    cudaFree(d_qmap);

    return qmap;
}


