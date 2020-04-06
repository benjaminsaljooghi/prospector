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


unsigned char* Prospector::get_qmap(ui* genome_encoding, ui genome_encoding_size)
{
    std::chrono::_V2::system_clock::time_point start;

    assert(K_START >= SIZE);

    cudaError er;

    ui bytes_genome_encoding = sizeof(ui) * genome_encoding_size;
    ui bytes_qmap = sizeof(uc) * genome_encoding_size * MAP_SIZE;
    
    ui* d_genome_encoding;
    uc* d_qmap;

    start = time();
    
    er = cudaMalloc(&d_genome_encoding, bytes_genome_encoding); checkCuda(er);
    er = cudaMemcpy(d_genome_encoding, genome_encoding, bytes_genome_encoding, cudaMemcpyHostToDevice); checkCuda(er);
    start = time(start, "genome encoding cudainit");

    er = cudaMalloc(&d_qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap malloc");

    er = cudaMemset(d_qmap, 0, bytes_qmap); checkCuda(er);
    start = time(start, "qmap memcpy");

    compute_qmap KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_genome_encoding, genome_encoding_size, d_qmap);

    uc* qmap;
    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap mallochost");

    cudaDeviceSynchronize();
    start = time(start, "kernel");

    er = cudaMemcpy(qmap, d_qmap, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er);
    start = time(start, "qmap memcpy");

    cudaFree(d_genome_encoding); cudaFree(d_qmap);

    return qmap;
}


