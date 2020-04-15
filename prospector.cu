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

#define DEBUG 1




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

void cudaWait()
{
    checkCuda ( cudaDeviceSynchronize() );
}

void Prospector::device_init()
{
    std::chrono::_V2::system_clock::time_point start = time();
    cudaFree(0);
    time(start, "device init");
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

Prospector::Encoding Prospector::get_genome_encoding(const char* genome, ui genome_size)
{
    auto start = time();

    cudaError er;
    
    char* d_genome;
    ui bytes_genome = sizeof(char) * genome_size;
    er = cudaMalloc(&d_genome, bytes_genome); checkCuda(er);
    er = cudaMemcpy(d_genome, genome, bytes_genome, cudaMemcpyHostToDevice); checkCuda(er);

    ui* encoding_d;
    ui encoding_size = genome_size - SIZE + 1;
    ui bytes_encoding = sizeof(ui) * encoding_size; 
    er = cudaMalloc(&encoding_d, bytes_encoding); checkCuda(er);

    compute_encoding KERNEL_ARGS3(128, 1024, 0) (d_genome, encoding_d, genome_size, encoding_size);

    ui* encoding;
    er = cudaMallocHost(&encoding, bytes_encoding); checkCuda(er);

    cudaWait();

    er = cudaMemcpy(encoding, encoding_d, bytes_encoding, cudaMemcpyDeviceToHost); checkCuda(er);

    time(start, "genome encoding total");

    return Encoding {
        encoding,
        encoding_d,
        encoding_size
    };
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

__global__ void compute_qmap_small(const ui* genome_encoding, const ui genome_encoding_size, uc* qmap)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    for (ui query = thread_id; query < genome_encoding_size - 200; query += stride) // TODO: replace 200 with MAP_SIZE or parameter map_size
    {
        ui query_enc = genome_encoding[query];
        for (ui i = 0; i < MAP_SIZE; i++)
        {
            ui t = genome_encoding[query + K_START + SPACER_SKIP + i];
            qmap[(query*MAP_SIZE) + i] = difference_gpu(query_enc, t);
        }
    }
}

__global__ void compute_qmap_big(const ui* encoding, const ui encoding_size, const ui* queries, const ui queries_size, uc* qmap, const ui map_size)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    for (ui q_i = thread_id; q_i < queries_size; q_i += stride) // warning. May go out of bounds if there is a query less than map_size away from the end of the genome
    {
        const ui query = queries[q_i];
        const ui query_enc = encoding[query];
        for (ui t_i = 0; t_i < map_size; t_i++)
        {
            // const ui t_enc = encoding[query + t_i]; // TODO: replace query + t_i with standard offset: query + K_START + SPACER_SKIP + i
            const ui t_enc = encoding[query + K_START + SPACER_SKIP + t_i]; // TODO: replace query + t_i with standard offset: query + K_START + SPACER_SKIP + i
            qmap[(q_i*map_size) + t_i] = difference_gpu(query_enc, t_enc);
        }
    }
}


uc* Prospector::get_qmap_small(ui* encoding_d, ui encoding_size)
{
    assert(K_START >= SIZE);
    cudaError er; 
    auto start = time();

    ui bytes_qmap = sizeof(uc) * encoding_size * MAP_SIZE;
    
    uc* qmap;
    uc* qmap_d;

    er = cudaMalloc(&qmap_d, bytes_qmap); checkCuda(er);
    start = time(start, "qmap small malloc");

    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap small mallochost");

    compute_qmap_small KERNEL_ARGS3(128, 1024, 0) (encoding_d, encoding_size, qmap_d);

    cudaWait();
    start = time(start, "qmap small kernel");

    er = cudaMemcpy(qmap, qmap_d, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er); 
    start = time(start, "qmap small memcpy");

    cudaFree(qmap_d);
    return qmap;
}


uc* Prospector::get_qmap_big(const ui* encoding_d, const ui encoding_size, const ui* queries, const ui queries_size, ui map_size)
{
    auto start = time();

    cudaError er; 

    ui bytes_queries = sizeof(ui) * queries_size;
    ui bytes_qmap = sizeof(uc) * queries_size * map_size;

    ui* queries_d;
    er = cudaMalloc(&queries_d, bytes_queries); checkCuda(er);
    start = time(start, "qmap big queries malloc");

    er = cudaMemcpy(queries_d, queries, bytes_queries, cudaMemcpyHostToDevice); checkCuda(er);
    start = time(start, "qmap big queries memcpy");

    uc* qmap;
    uc* qmap_d;

    er = cudaMalloc(&qmap_d, bytes_qmap); checkCuda(er);
    start = time(start, "qmap big malloc");

    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
    start = time(start, "qmap big mallochost");

    compute_qmap_big KERNEL_ARGS3(128, 1024, 0) (encoding_d, encoding_size, queries_d, queries_size, qmap_d, map_size);
    
    cudaWait();
    start = time(start, "qmap big kernel");

    er = cudaMemcpy(qmap, qmap_d, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er); 
    start = time(start, "qmap big memcpy");

    cudaFree(qmap_d);
    return qmap;
}


