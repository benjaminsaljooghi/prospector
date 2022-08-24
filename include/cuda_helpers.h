#pragma once

#ifndef PROSPECTOR_CUDA_H
#define PROSPECTOR_CUDA_H

#include "cuda.h"

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
#define GRID 32
#define BLOCK 256

cudaError_t checkCudaAlways(cudaError_t result);
cudaError_t checkCuda(cudaError_t result);
void cudaWait();

#endif //PROSPECTOR_CUDA_H
