#include "prospector.h"
#include "cuda_helpers.h"
#include <cassert>
#include <cstdio>

__device__ int __popc (ui x);
__device__ uc difference_gpu(const ui& _a, const ui& _b)
{
    ui _xor = (_a ^ _b);
    ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ui oddBits = _xor & 0x5555555555555555ull;
    ui comp = (evenBits >> 1) | oddBits;
    return __popc(comp);
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
    
    printf("have invalid char: %c\n", c);
    assert(0);
    //abort();
    return 0;
}

__global__ void compute_encoding(const char* genome, kmer* encoding, ull genome_size, ull encoding_size)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    for (ui i = thread_id; i < encoding_size; i += stride)
    {
        ui e = 0;
        for (int j = 0; j < Prospector::size; j++)
        {
            e |= scheme(genome[i + j]) << (j * Prospector::bits);
        }
        encoding[i] = e;
    }
}

Prospector::Encoding Prospector::get_genome_encoding(const char* genome, ull genome_size)
{

    Prospector::Encoding encoding{};

    cudaError er;
    
    char* d_genome;
    ull bytes_genome = sizeof(char) * genome_size;
    er = cudaMalloc(&d_genome, bytes_genome); checkCuda(er);
    er = cudaMemcpy(d_genome, genome, bytes_genome, cudaMemcpyHostToDevice); checkCuda(er);

    encoding.size = genome_size - Prospector::size + 1;
    encoding.bytes = sizeof(ui) * encoding.size;
    er = cudaMalloc(&encoding.d, encoding.bytes); checkCuda(er);

    compute_encoding KERNEL_ARGS3(GRID, BLOCK, 0) (d_genome, encoding.d, genome_size, encoding.size);

    er = cudaMallocHost(&encoding.h, encoding.bytes); checkCuda(er);

    cudaWait();

    er = cudaMemcpy(encoding.h, encoding.d, encoding.bytes, cudaMemcpyDeviceToHost); checkCuda(er);

    cudaFree(d_genome);

    return encoding;
}

__global__ void compute_qmap_small(const ui* encoding, const ui encoding_size, uc* qmap)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    for (ui query = thread_id; query < encoding_size - 200; query += stride) // TODO: replace 200 with MAP_SIZE or parameter map_size
    {
        ui query_enc = encoding[query];
        for (ui i = 0; i < Prospector::map_size_small; i++)
        {
            ui t = encoding[query + Prospector::k_start + Prospector::spacer_skip + i];
            qmap[(query * Prospector::map_size_small) + i] = difference_gpu(query_enc, t);
        }
    }
}

//__global__ void compute_qmap_big(const ui* encoding, const ui encoding_size, const ui* queries, const ui queries_size, uc* qmap)
//{
//    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
//    const ui stride = blockDim.x * gridDim.x;
//
//    for (ui q_i = thread_id; q_i < queries_size; q_i += stride) // warning. May go out of bounds if there is a query less than map_size away from the end of the genome
//    {
//        const ui query = queries[q_i];
//        const ui query_enc = encoding[query];
//        for (ui t_i = 0; t_i < Prospector::map_size_big; t_i++)
//        {
//            const ui t_enc = encoding[query + Prospector::k_start + Prospector::spacer_skip + t_i];
//            qmap[(q_i * Prospector::map_size_big) + t_i] = difference_gpu(query_enc, t_enc);
//        }
//    }
//}



uc* qmap;
uc* qmap_d;

void Prospector::device_init()
{
    cudaDeviceReset();
    cudaError_t error = checkCudaAlways(cudaFree(0));
    if (error != cudaSuccess)
    {
        printf("device initialization error\n");
        abort();
    }

    cudaError er;
    er = cudaMalloc(&qmap_d, 2000000000 ); checkCuda(er);

    er = cudaMallocHost(&qmap, 2000000000); checkCuda(er);

}

uc* Prospector::get_qmap_small(const kmer* encoding, const ull encoding_size)
{
    assert(Prospector::k_start >= Prospector::size);
    cudaError er; 

    ui bytes_qmap = sizeof(uc) * encoding_size * Prospector::map_size_small;

    printf("bytes_qmap: %d\n", bytes_qmap);
    
    compute_qmap_small KERNEL_ARGS3(GRID, BLOCK, 0) (encoding, encoding_size, qmap_d);

    cudaWait();

    er = cudaMemcpy(qmap, qmap_d, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er); 

    //cudaFree(qmap_d);
    return qmap;
}


//uc* Prospector::get_qmap_big(const ui* encoding_d, const ui encoding_size, const ui* queries, const ui queries_size)
//{
//    auto start = time();
//
//    cudaError er; 
//
//    ui bytes_queries = sizeof(ui) * queries_size;
//    ui bytes_qmap = sizeof(uc) * queries_size * Prospector::map_size_big;
//
//    ui* queries_d;
//    er = cudaMalloc(&queries_d, bytes_queries); checkCuda(er);
//    start = time(start, "qmap big queries malloc");
//
//    er = cudaMemcpy(queries_d, queries, bytes_queries, cudaMemcpyHostToDevice); checkCuda(er);
//    start = time(start, "qmap big queries memcpy");
//
//    uc* qmap;
//    uc* qmap_d;
//
//    er = cudaMalloc(&qmap_d, bytes_qmap); checkCuda(er);
//    start = time(start, "qmap big malloc");
//
//    er = cudaMallocHost(&qmap, bytes_qmap); checkCuda(er);
//    start = time(start, "qmap big mallochost");
//
//    compute_qmap_big KERNEL_ARGS3(GRID, BLOCK, 0) (encoding_d, encoding_size, queries_d, queries_size, qmap_d);
//    
//    cudaWait();
//    start = time(start, "qmap big kernel");
//
//    er = cudaMemcpy(qmap, qmap_d, bytes_qmap, cudaMemcpyDeviceToHost); checkCuda(er); 
//    start = time(start, "qmap big memcpy");
//
//    cudaFree(qmap_d);
//    return qmap;
//}

void Prospector::free_encoding(Encoding encoding)
{
    cudaFree(encoding.d);
    cudaFreeHost(encoding.h);
}



