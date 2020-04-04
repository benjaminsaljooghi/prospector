#include "prospector.h"


#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"


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










typedef unsigned long long ull;

#define BITS 4
#define SIZE 16
map<char, ull> scheme {
    {'A', 1},
    {'C', 2},
    {'G', 4},
    {'T', 8} 
};


ull encoded(const string& kmer)
{
    assert(kmer.size() == SIZE);
    ull e = 0;
    for (int i = 0; i < kmer.size(); i++)
        e += scheme.at(kmer[i]) << (i * BITS);
    return e;
}

vector<ull> encoded_genome(const string& genome)
{
    double __start = omp_get_wtime();
    vector<ull> encoding;
    for (size_t i = 0; i < genome.size() - SIZE + 1; i++)
    {
        string kmer = genome.substr(i, SIZE);
        encoding.push_back(encoded(kmer));
    }
    done(__start, "\t\tencoding");
    return encoding;
}

__device__ ull difference(const ull& _a, const ull& _b)
{
    ull _xor = _a ^ _b;
    ull _and = _xor & _a;
    ull _pc = __popcll(_and);
    // printf("%d\n", _pc);
    return _pc;
}

__device__ bool mutant(const char* genome, const ull* genome_encoding, const unsigned int& i, const unsigned int& j, const unsigned int& k, const unsigned int& allowed_mutations)
{
    ull diff = 0;
    const int chunks = k / SIZE;
    for (int chunk = 0; chunk < chunks; chunk++)
    {
        diff += difference(genome_encoding[i + (chunk * SIZE)], genome_encoding[j + (chunk * SIZE)]);
        if (diff > allowed_mutations)
        {
            return false;
        }
    }

    return true;

    // compute final diffs
    const int chars = k - (chunks * SIZE);
    for (int z = 0; z < chars; z++)
    {
        diff += genome[i + (chunks * SIZE) +  z] == genome[j + (chunks * SIZE) +  z] ? 0 : 1;
    }
    return diff <= allowed_mutations;
    return true;
}









void cwait()
{
	cudaError err = cudaDeviceSynchronize();
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(err);
	}
}



inline
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




#define C_GRID 64
#define C_BLOCK 512
#define THREAD_BUFFER_COUNT 150 // this can be brought lower and lower until FATAL messages are received

// #define C_GRID 1
// #define C_BLOCK 1
// #define THREAD_BUFFER_COUNT (150 * 64 * 512) 

#define K_BUFFER_COUNT (THREAD_BUFFER_COUNT * C_GRID * C_BLOCK)
#define TOTAL_BUFFER_COUNT (K_COUNT * K_BUFFER_COUNT)


__global__ void discover_crisprs(
        const char* genome,
        const ull* genome_encoding,
        const ull genome_encoding_size,
        unsigned int* buffer,
        unsigned int* starts,
        unsigned char* sizes,
        unsigned int buffer_start,
        unsigned int query_start,
        unsigned int i
    )
{
    const ull thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ull stride = blockDim.x * gridDim.x;
    const ull k = K_START + i;
    const ull allowed_mutations = k / MUTANT_TOLERANCE_RATIO;
    const ull thread_buffer_start = buffer_start + (thread_id * THREAD_BUFFER_COUNT);
    const ull thread_buffer_limit = thread_buffer_start + THREAD_BUFFER_COUNT;

    ull buffer_pointer = thread_buffer_start;
    for (ull query = thread_id; query < genome_encoding_size; query += stride)
    {
        ull bound = query + k + SPACER_SKIP;
        starts[query] = buffer_pointer;
        buffer[buffer_pointer++] = query;
        for (ull target = bound; target - bound <= SPACER_MAX && target < genome_encoding_size; target++)
        {
            if (mutant(genome, genome_encoding, query, target, k, allowed_mutations))
            {
                buffer[buffer_pointer++] = target;
                bound = target + k + SPACER_SKIP;
                target = bound;
            }
        }

        sizes[query] = buffer_pointer - starts[query];

        // ull size = sizes[query];
        // if (size>5)
        // {
        //     printf("%lu\n", size);
        //     for (int x = 0; x < size; x++)
        //     {
        //         ull start = starts[query];
        //         printf("%lu %lu\n", start + x, buffer[start + x]);
        //     }
        //     printf("\n");
        // }



        #if DEBUG == 1
            if (buffer_pointer  >= thread_buffer_limit)
            {
                printf("FATAL: exceeded thread buffer limit. %lu/%lu\n", buffer_pointer, thread_buffer_limit);
            }
        #endif

    }

}

vector<Crispr> crispr_formation_single(const char* device_genome, const ull* device_genome_encoding, const ull genome_size, const ull genome_encoding_size, const ull i)
{
    ull k = K_START + i;

    cudaError er;

    ull count_buffer = K_BUFFER_COUNT;
    ull count_starts = genome_size - k + 1; // will form (at least) a singleton crispr for each i in the genome
    ull count_sizes = genome_size - k + 1;

    ull bytes_buffer = 4 * count_buffer;
    ull bytes_starts = 4 * count_starts;
    ull bytes_sizes = 1 * count_sizes;

    // buffer info
    unsigned int* buffer, *device_buffer;
    unsigned int* starts, *device_starts;
    unsigned char* sizes, *device_sizes;

    er = cudaMallocHost(&buffer, bytes_buffer); checkCuda(er);
    er = cudaMallocHost(&starts, bytes_starts); checkCuda(er);
    er = cudaMallocHost(&sizes, bytes_sizes); checkCuda(er);

    er = cudaMalloc(&device_buffer, bytes_buffer); checkCuda(er);
    er = cudaMalloc(&device_starts, bytes_starts); checkCuda(er);
    er = cudaMalloc(&device_sizes, bytes_sizes); checkCuda(er);

    er = cudaMemset(device_buffer, 0, bytes_buffer); checkCuda(er);
    er = cudaMemset(device_starts, 0, bytes_starts); checkCuda(er);
    er = cudaMemset(device_sizes, 0, bytes_sizes); checkCuda(er);

    printf("\t\tkernel... ");
    double start = omp_get_wtime();
    discover_crisprs KERNEL_ARGS3(C_GRID, C_BLOCK, 0) 
    (
        device_genome,
        device_genome_encoding,
        genome_encoding_size,
        device_buffer,
        device_starts,
        device_sizes,
        0,
        0,
        i
    );

    cwait();

    done(start);

    er = cudaMemcpy(buffer, device_buffer, bytes_buffer, cudaMemcpyDeviceToHost); checkCuda(er);
    er = cudaMemcpy(starts, device_starts, bytes_starts, cudaMemcpyDeviceToHost); checkCuda(er);
    er = cudaMemcpy(sizes, device_sizes, bytes_sizes, cudaMemcpyDeviceToHost); checkCuda(er);


    start = omp_get_wtime();
    
    printf("\t\textract...");
    vector<Crispr> crisprs;
    for (ull q = 0; q < count_starts; q++)
    {
        unsigned char __size = *(sizes + q);
        if (__size >= MIN_REPEATS)
        {
            unsigned int start = *(starts + q);
            vector<unsigned int> genome_indices;
            for (int walk = 0; walk < __size; walk++)
            {
                genome_indices.push_back(  *(buffer + start + walk)  ); // just use indices/ptrs rather than a vector?
            }
            assert(__size == genome_indices.size()); // not necessary, just put it in to remind you later of double up on size in constructor
            Crispr crispr(k, genome_indices, genome_indices.size()); // why the double up here..? refactor crispr?
            crisprs.push_back(crispr);
        }
    }
    done(start);

    return crisprs;

    // vector<Crispr> crisprs;
    // unsigned int __dyad_start = 0;
    // for (unsigned int i = 0; i < K_COUNT; i++)
    // {
    //     unsigned int k = K_START + i;
    //     vector<Crispr> crisprs;

    //     for (unsigned int dyad_index = 0; dyad_index < dyad_counts[i]; dyad_index++)
    //     {
    //         unsigned char __size = *(sizes + __dyad_start + dyad_index);
    //         if (__size >= MIN_REPEATS)
    //         {
    //             unsigned int start = *(starts + __dyad_start + dyad_index);
    //             unsigned int* _s = buffer + start;
    //             unsigned int* _e = _s + __size;
    //             Crispr crispr(k, _s, _e);
    //             crisprs.push_back(crispr);
    //         }
    //     }

    //     crisprs.insert(crisprs.end(), crisprs.begin(), crisprs.end());
    //     __dyad_start += dyad_counts[i];
    // }


}

vector<Crispr> crispr_formation(const char* device_genome, const ull* device_genome_encoding, const ull genome_size, const ull genome_encoding_size)
{

    vector<Crispr> all_crisprs;
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        vector<Crispr> crisprs = crispr_formation_single(device_genome, device_genome_encoding, genome_size, genome_encoding_size, i);
        for (ull j = 0; j < crisprs.size(); j++)
        {
            all_crisprs.push_back(crisprs[j]);
        }
    }

    return all_crisprs;
}



vector<Crispr> prospector_main_gpu(const string& genome)
{

    cudaDeviceReset();

    vector<ull> genome_encoding = encoded_genome(genome);

    char* device_genome;
    ull* device_genome_encoding;

    cudaMalloc(&device_genome, 1 * genome.length());
    cudaMalloc(&device_genome_encoding, 8 * genome_encoding.size());

    cudaMemcpy(device_genome, genome.c_str(), 1 * genome.length(), cudaMemcpyHostToDevice);
    cudaMemcpy(device_genome_encoding, &genome_encoding[0], 8 * genome_encoding.size(), cudaMemcpyHostToHost);

    // start = omp_get_wtime();
     // printf("\tdyad discovery\n");
    // Dyad_Info dyad_info = dyad_discovery(d_genome, genome.length());
    // done(start, "\tdyad discovery");

    printf("\tcrispr formation\n");
    double start = omp_get_wtime();
        vector<Crispr> crisprs = crispr_formation(device_genome, device_genome_encoding, genome.length(), genome_encoding.size());
    done(start, "\tcrispr formation");

    cudaFree(device_genome);

    return crisprs;
}


vector<Crispr> Prospector::prospector_main(const string& genome)
{
    printf("genome has size %zd\n", genome.size());
    
    double start;
    printf("prospector\n");
    start = omp_get_wtime();
    vector<Crispr> crisprs = prospector_main_gpu(genome);
    done(start, "prospector");

    return crisprs;
}


