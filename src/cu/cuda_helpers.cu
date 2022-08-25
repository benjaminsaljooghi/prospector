#include "cuda_helpers.h"
#include <cstdio>

cudaError_t checkCudaAlways(cudaError_t result)
{
    if (result != cudaSuccess)
    {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    }
    return result;
}

cudaError_t checkCuda(cudaError_t result)
{
#if DEBUG == 1
    checkCudaAlways(result);
#endif
    return result;
}

void cudaWait()
{
    checkCuda ( cudaDeviceSynchronize() );
}
