#include "prospector.h"


#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"

#include <cassert>



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


unsigned char* Prospector::get_qmap(ui* genome_encoding, ui genome_encoding_size)
{
    assert(K_START >= SIZE);

    cudaError er;

    ui* d_genome_encoding;

    er = cudaMalloc(&d_genome_encoding, 4 * genome_encoding_size); checkCuda(er);
    er = cudaMemcpy(d_genome_encoding, &genome_encoding[0], 4 * genome_encoding_size, cudaMemcpyHostToDevice); checkCuda(er);
    ui count_qmap = genome_encoding_size * MAP_SIZE;
    ui bytes_qmap = 1 * count_qmap;
    unsigned char* qmap, *d_qmap;
    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
    er = cudaMalloc(&d_qmap, bytes_qmap); checkCuda(er);
    er = cudaMemset(d_qmap, 0, bytes_qmap); checkCuda(er);

    compute_qmap KERNEL_ARGS3(C_GRID, C_BLOCK, 0) (d_genome_encoding, genome_encoding_size, d_qmap);

    cudaDeviceSynchronize();

    er = cudaMemcpy(qmap, d_qmap, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er);

    cudaFree(d_genome_encoding); cudaFree(d_qmap);

    return qmap;
}


