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


void cwait()
{
	cudaError err = cudaDeviceSynchronize();
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(err);
	}
}

__device__ char complement(char nuc)
{
    switch (nuc)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    case 'N':
        return 'N';
    case 'n':
        return 'n';
    default:
        return 'n';
    }
}

__device__ bool mutant(const char* genome, unsigned int start_a, unsigned int start_b, unsigned int k)
{
	unsigned int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;

	unsigned int mutations = 0;

	for (unsigned int i = 0; i < k; i++)
	{
        mutations += genome[start_a + i] == genome[start_b + i] ? 0 : 1;
		if (mutations > allowed_mutations)
        {
			return false;
        }
	}
	return true;
}

#if DEBUG == 1
__device__ bool is_dyad_debug_check(unsigned int start_index)
{
    return start_index >= DEBUG_START && start_index <= DEBUG_END;
}
#endif


__device__ bool is_dyad(const char* genome, unsigned int start_index, unsigned int k)
{
#if DEBUG == 1
    if (!is_dyad_debug_check(start_index))
    {
        return false;
    }
#endif

    unsigned int end_index = start_index + k - 1;

    unsigned int range = k/2 - 1;
    unsigned int mismatch_count = 0;
    for (unsigned int i = 0; i < range; i++)
	{
		char upstream = genome[start_index + i];
		char downstream = genome[end_index - i];
        mismatch_count += upstream == complement(downstream) ? 0 : 1;
	}

    double mismatch_ratio = (double) mismatch_count / (double) range;
    return mismatch_ratio < 0.8;
}



#define C_GRID 64
#define C_BLOCK 512
#define THREAD_BUFFER_COUNT 100 // this can be brought lower and lower until FATAL messages are received
#define K_BUFFER_COUNT (THREAD_BUFFER_COUNT * C_GRID * C_BLOCK)
#define TOTAL_BUFFER_COUNT (K_COUNT * K_BUFFER_COUNT)


__global__ void discover_crisprs(const char* genome, unsigned int* dyads, const size_t* dyad_counts, unsigned int* buffer, unsigned int* starts, unsigned char* sizes, unsigned int buffer_start, unsigned int dyad_start, unsigned int i)
{
    const unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = blockDim.x * gridDim.x;
    const unsigned int dyad_count = dyad_counts[i];
    const unsigned int k = K_START + i;
    const unsigned int thread_buffer_start = buffer_start + (thread_id * THREAD_BUFFER_COUNT);
    const unsigned int thread_buffer_limit = thread_buffer_start + THREAD_BUFFER_COUNT;

    unsigned int buffer_pointer = thread_buffer_start;
    for (unsigned int query_d_index = thread_id; query_d_index < dyad_count; query_d_index += stride)
    {
        unsigned int query_dyad = *(dyads + dyad_start + query_d_index);
        unsigned int bound = query_dyad + k + SPACER_SKIP;

        *(starts + dyad_start + query_d_index) = buffer_pointer;
        *(buffer + buffer_pointer) = query_dyad;
        unsigned char crispr_size = 1;

        for (unsigned int target_d_index = query_d_index + 1; target_d_index < dyad_count; target_d_index++) // this for loop goes up to the dyad count but in practice this will never happen. May want to rewrite this. while loop? or for loop up to CRISPR_BUFFER?
        {
            unsigned int target_dyad = *(dyads + dyad_start + target_d_index);
            if (target_dyad < bound) continue;
            if (target_dyad - bound > SPACER_MAX) break;
            if (mutant(genome, query_dyad, target_dyad, k))
            {
                *(buffer + buffer_pointer + crispr_size++) = target_dyad;
                bound = target_dyad + k + SPACER_SKIP;
            }
        }

        buffer_pointer += crispr_size;
        *(sizes + dyad_start + query_d_index) = crispr_size;

#if DEBUG == 1
        if (buffer_pointer  >= thread_buffer_limit)
        {
            printf("FATAL: exceeded thread buffer limit. %d/%d\n", buffer_pointer, thread_buffer_limit);
        }
#endif
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


struct Dyad_Info
{
    vector<unsigned int> dyads;
    size_t* dyad_counts;
    size_t total_dyad_count;
};


vector<Crispr> crispr_formation(const char* device_genome, Dyad_Info dyad_info)
{
    double start_time;

    unsigned int* dyads = &dyad_info.dyads[0];
    size_t* dyad_counts = dyad_info.dyad_counts;
    size_t total_dyad_count = dyad_info.total_dyad_count;

    cudaError er;

    size_t count_dyads = total_dyad_count;
    size_t count_dyad_counts = K_COUNT;
    size_t count_buffer = TOTAL_BUFFER_COUNT;
    size_t count_starts = total_dyad_count;
    size_t count_sizes = total_dyad_count;

    size_t bytes_dyads = 4 * count_dyads;
    size_t bytes_dyad_counts = 8 * count_dyad_counts;
    size_t bytes_buffer = 4 * count_buffer;
    size_t bytes_starts = 4 * count_starts;
    size_t bytes_sizes = 1 * count_sizes;

    // dyad info
    unsigned int* d_dyads;
    size_t* d_dyad_counts;
    er = cudaMalloc(&d_dyads, bytes_dyads); checkCuda(er);
    er = cudaMalloc(&d_dyad_counts, bytes_dyad_counts); checkCuda(er);
    start_time = omp_get_wtime();
    er = cudaMemcpy(d_dyads, dyads, bytes_dyads, cudaMemcpyHostToDevice); checkCuda(er);
    er = cudaMemcpy(d_dyad_counts, dyad_counts, bytes_dyad_counts, cudaMemcpyHostToDevice); checkCuda(er);
    done(start_time, "\t\tmemcpies");

    // buffer info
    unsigned int* buffer, *d_buffer;
    unsigned int* starts, *d_starts;
    unsigned char* sizes, *d_sizes;

    er = cudaMallocHost(&buffer, bytes_buffer); checkCuda(er);
    er = cudaMallocHost(&starts, bytes_starts); checkCuda(er);
    er = cudaMallocHost(&sizes, bytes_sizes); checkCuda(er);

    er = cudaMalloc(&d_buffer, bytes_buffer); checkCuda(er);
    er = cudaMalloc(&d_starts, bytes_starts); checkCuda(er);
    er = cudaMalloc(&d_sizes, bytes_sizes); checkCuda(er);

    start_time = omp_get_wtime();
    er = cudaMemset(buffer, 0, bytes_buffer); checkCuda(er);
    er = cudaMemset(starts, 0, bytes_starts); checkCuda(er);
    er = cudaMemset(sizes, 0, bytes_sizes); checkCuda(er);
    done(start_time, "\t\tmemsets");

    // stream init
    unsigned int buffer_starts[K_COUNT];
    unsigned int dyad_starts[K_COUNT];
    buffer_starts[0] = 0; dyad_starts[0] = 0;
    for (int i = 1; i < K_COUNT; i++)
    {
        buffer_starts[i] = buffer_starts[i-1] + K_BUFFER_COUNT;
        dyad_starts[i] = dyad_starts[i-1] + dyad_counts[i-1];
    }

    cudaStream_t stream[K_COUNT];
    for (unsigned int i = 0; i < K_COUNT; i++) cudaStreamCreate(&stream[i]);

    start_time = omp_get_wtime();
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        discover_crisprs KERNEL_ARGS4(C_GRID, C_BLOCK, 0, stream[i]) (device_genome, d_dyads, d_dyad_counts, d_buffer, d_starts, d_sizes, buffer_starts[i], dyad_starts[i], i);
        er = cudaMemcpyAsync(&buffer[buffer_starts[i]], &d_buffer[buffer_starts[i]], K_BUFFER_COUNT * 4, cudaMemcpyDeviceToHost, stream[i]); checkCuda(er);
        er = cudaMemcpyAsync(&starts[dyad_starts[i]], &d_starts[dyad_starts[i]], dyad_counts[i] * 4, cudaMemcpyDeviceToHost, stream[i]); checkCuda(er);
        er = cudaMemcpyAsync(&sizes[dyad_starts[i]], &d_sizes[dyad_starts[i]], dyad_counts[i] * 1, cudaMemcpyDeviceToHost, stream[i]); checkCuda(er);
    }

    cwait();
    done(start_time, "\t\tasync kernels and memcpies");

    start_time = omp_get_wtime();
    vector<Crispr> all_crisprs;
    unsigned int __dyad_start = 0;
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        unsigned int k = K_START + i;
        vector<Crispr> crisprs;

        for (unsigned int dyad_index = 0; dyad_index < dyad_counts[i]; dyad_index++)
        {
            unsigned char __size = *(sizes + __dyad_start + dyad_index);
            if (__size >= MIN_REPEATS)
            {
                unsigned int start = *(starts + __dyad_start + dyad_index);
                unsigned int* _s = buffer + start;
                unsigned int* _e = _s + __size;
                Crispr crispr(k, _s, _e);
                crisprs.push_back(crispr);
            }
        }

        all_crisprs.insert(all_crisprs.end(), crisprs.begin(), crisprs.end());
        __dyad_start += dyad_counts[i];
    }

    done(start_time, "\t\textract");


    return all_crisprs;
}


__global__ void dyad_discovery(const char* genome, size_t genome_len, bool* dyad_buffer)
{
    const unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = blockDim.x * gridDim.x;

    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        const unsigned int k = K_START + i;
        bool* buffer_pointer = dyad_buffer + (i * genome_len);
        for (unsigned int dyad = thread_id; dyad + k < genome_len; dyad += stride)
        {
            *(buffer_pointer + dyad) = is_dyad(genome, dyad, k);
        }
    }
}


Dyad_Info dyad_discovery(char* device_genome, size_t genome_len)
{
    double start_time;

    size_t buffer_count = genome_len * K_COUNT;
    size_t buffer_bytes = 1 * buffer_count;

    bool *buffer, *d_buffer;
    cudaMallocHost(&buffer, buffer_bytes);
    cudaMalloc(&d_buffer, buffer_bytes);
    cudaMemset(d_buffer, 0, buffer_bytes);

    start_time = omp_get_wtime();
    dyad_discovery KERNEL_ARGS2(128, 512) (device_genome, genome_len, d_buffer);
    cwait();
    done(start_time, "\t\tkernel");

    start_time = omp_get_wtime();
    cudaMemcpy(buffer, d_buffer, buffer_bytes, cudaMemcpyDeviceToHost);
    cudaFree(d_buffer);
    done(start_time, "\t\tmemcpy");

    start_time = omp_get_wtime();
    vector<unsigned int> dyads;
    size_t* dyad_counts = (size_t*) malloc(K_COUNT * 8);
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        bool* buffer_pointer = buffer + (i * genome_len);
        size_t count_before = dyads.size();
        for (unsigned int j = 0; j < genome_len; j++)
        {
            if (*(buffer_pointer++) == 1)
            {
                dyads.push_back(j);
            }
        }
        dyad_counts[i] = dyads.size() - count_before;
    }
    done(start_time, "\t\textract");

    cudaFreeHost(buffer);
    cudaFree(d_buffer);

    Dyad_Info info;
    info.dyads = dyads;
    info.dyad_counts = dyad_counts;
    info.total_dyad_count = dyads.size();
    return info;
}


vector<Crispr> prospector_main_gpu(const string& genome)
{
    cudaDeviceReset();

    double start;

    char* d_genome;
    cudaMalloc(&d_genome, 1 * genome.length());
    cudaMemcpy(d_genome, genome.c_str(), 1 * genome.length(), cudaMemcpyHostToDevice);

    start = omp_get_wtime();
    printf("\tdyad discovery\n");
    Dyad_Info dyad_info = dyad_discovery(d_genome, genome.length());
    done(start, "\tdyad discovery");


    printf("\tcrispr formation\n");
    start = omp_get_wtime();
    vector<Crispr> crisprs = crispr_formation(d_genome, dyad_info);
    done(start, "\tcrispr formation");

    cudaFree(d_genome);

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


