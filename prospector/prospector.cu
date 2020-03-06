// CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"


//#include "../util/stdafx.h"
//#include "../util/util.h"
#include "prospector.h"


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
	printf("waiting for kernel... ");
	clock_t start = clock();
	cudaError err = cudaDeviceSynchronize();
	printf("done in %.3f seconds\n", duration(start));
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(err);
	}
}

void cufree(void* device_ptr)
{
	printf("executing cudafree... ");
    clock_t start = clock();
	cudaError err = cudaFree(device_ptr);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to free device ptr (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);
}

template <typename T> void cpull(T* host, const T* device, unsigned int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;

	printf("memcpy %*zd bytes from device... ", printf_BYTE_FORMAT_ALIGN, bytes);
    clock_t start = clock();
	err = cudaMemcpy(host, device, bytes, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from device to host (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);
    
}

template <typename T> T* cpush(const T* src, unsigned int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;
	T* ptr = nullptr;

	printf("malloc+memcpy %*zd bytes to device... ", printf_BYTE_FORMAT_ALIGN, bytes);
    clock_t start = clock();
	err = cudaMalloc((void**)& ptr, bytes);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to malloc device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    // done(start);

	// printf("memcpy %*zd bytes to device... ", printf_BYTE_FORMAT_ALIGN, bytes);
    // start = clock();
	err = cudaMemcpy(ptr, src, bytes, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);

	return (T*)ptr;
}

__device__ char complement(char nuc)
{
    // printf("%c\n", nuc);
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

__device__ bool is_dyad_debug_check(unsigned int start_index)
{
    return start_index >= DEBUG_START && start_index <= DEBUG_END;
}

__device__ bool is_dyad(const char* genome, unsigned int start_index, unsigned int k)
{
    if (!is_dyad_debug_check(start_index))
    {
        return false;
    }

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
    return mismatch_ratio < 0.75;
}


int atomicAdd(int* address, int val);
unsigned int atomicAdd(unsigned int* address,
                       unsigned int val);
unsigned long long int atomicAdd(unsigned long long int* address,
                                 unsigned long long int val);
float atomicAdd(float* address, float val);
double atomicAdd(double* address, double val);
__half2 atomicAdd(__half2 *address, __half2 val);
__half atomicAdd(__half *address, __half val);



#define C_GRID 16
#define C_BLOCK 256
#define THREAD_BUFFER_COUNT 1000
#define K_BUFFER_COUNT (THREAD_BUFFER_COUNT * C_GRID * C_BLOCK)
#define TOTAL_BUFFER_COUNT (K_COUNT * K_BUFFER_COUNT)

__global__ void discover_crisprs(const char* genome, unsigned int* dyads, const size_t* dyad_sizes, unsigned int* buffer, unsigned int* starts, unsigned char* sizes)
{
    unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;


    unsigned int dyad_start = 0;

    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        unsigned int k = K_START + i;

        const unsigned int k_buffer_start = i * K_BUFFER_COUNT;
        const unsigned int thread_buffer_start = k_buffer_start + (thread_id * THREAD_BUFFER_COUNT);

        const unsigned int thread_buffer_limit = thread_buffer_start + THREAD_BUFFER_COUNT; // exclusive limit

        unsigned int buffer_pointer = thread_buffer_start;

        for (unsigned int query_d_index = thread_id; query_d_index < dyad_sizes[i]; query_d_index += stride)
        {
            unsigned int query_dyad = *(dyads + dyad_start + query_d_index);
            unsigned int bound = query_dyad + k + SPACER_SKIP;

            *(starts + dyad_start + query_d_index) = buffer_pointer;
            *(buffer + buffer_pointer) = query_dyad;
            unsigned int crispr_size = 1;

            for (unsigned int target_d_index = query_d_index + 1; target_d_index < dyad_sizes[i]; target_d_index++) // this for loop goes up to the dyad count but in practice this will never happen. May want to rewrite this. while loop? or for loop up to CRISPR_BUFFER?
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

            if (buffer_pointer  >= thread_buffer_limit)
            {
                printf("FATAL: exceeded thread buffer limit. %d/%d\n", buffer_pointer, thread_buffer_limit);
//                exit(1);
            }
        }

        dyad_start += dyad_sizes[i];
    }


//    printf("final buffer_pointer: %d\n", (buffer_pointer - initial_buffer_pointer));
}



vector<Crispr> crispr_gen(const char* device_genome, vector<vector<unsigned int>> all_dyads)
{
    clock_t start_time;
    vector<Crispr> all_crisprs;

    auto* dyad_sizes = (size_t*) malloc(K_COUNT * sizeof(size_t));

    size_t total_dyad_count = 0;
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        total_dyad_count += all_dyads[i].size();
        dyad_sizes[i] = all_dyads[i].size();
    }

    vector<unsigned int> flat = flatten(all_dyads);

    size_t bytes_buffer = 4 * TOTAL_BUFFER_COUNT;
    size_t bytes_dyads = 4 * total_dyad_count;
    size_t bytes_starts = 4 * total_dyad_count;
    size_t bytes_sizes = 1 * total_dyad_count;

    auto* buffer = (unsigned int*) malloc(bytes_buffer);
    auto* dyads = &flat[0];
    auto* starts = (unsigned int*) malloc(bytes_starts);
    auto* sizes = (unsigned char*) malloc(bytes_sizes);

    unsigned int* d_buffer;
    size_t* d_dyad_sizes;
    unsigned int* d_dyads;
    unsigned int* d_starts;
    unsigned char* d_sizes;

    // cuda malloc
    auto erri = cudaMalloc(&d_buffer, bytes_buffer);
    auto errj = cudaMalloc(&d_dyads, bytes_dyads);
    auto errk = cudaMalloc(&d_dyad_sizes, K_COUNT * sizeof(size_t));
    auto errl = cudaMalloc(&d_starts, bytes_starts);
    auto errm = cudaMalloc(&d_sizes, bytes_sizes);

    // cuda memcpy
    auto erra = cudaMemcpy(d_buffer, buffer, bytes_buffer, cudaMemcpyHostToDevice);
    auto errb = cudaMemcpy(d_dyads, dyads, bytes_dyads, cudaMemcpyHostToDevice);
    auto errc = cudaMemcpy(d_dyad_sizes, dyad_sizes, K_COUNT * sizeof(size_t), cudaMemcpyHostToDevice);
    auto errd = cudaMemcpy(d_starts, starts, bytes_starts, cudaMemcpyHostToDevice);
    auto erre = cudaMemcpy(d_sizes, sizes, bytes_sizes, cudaMemcpyHostToDevice);


    discover_crisprs KERNEL_ARGS2(C_GRID, C_BLOCK) (device_genome, d_dyads, d_dyad_sizes, d_buffer, d_starts, d_sizes);
    cwait();

    cudaMemcpy(buffer, d_buffer, bytes_buffer, cudaMemcpyDeviceToHost);
    cudaMemcpy(starts, d_starts, bytes_starts, cudaMemcpyDeviceToHost);
    cudaMemcpy(sizes, d_sizes, bytes_sizes, cudaMemcpyDeviceToHost);


//    cufree(d_dyads);
//    cufree(d_dyad_sizes);
//    cufree(d_buffers);
//    cufree(d_sizes);

    printf("extract..."); start_time = clock();

    unsigned int dyad_start = 0;
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        unsigned int k = K_START + i;
        vector<Crispr> crisprs;

        for (unsigned int dyad_index = 0; dyad_index < dyad_sizes[i]; dyad_index++)
        {
            unsigned char __size = *(sizes + dyad_start + dyad_index);
            if (__size >= MIN_REPEATS)
            {
                unsigned int start = *(starts + dyad_start + dyad_index);
                unsigned int* _s = buffer + start;
                unsigned int* _e = _s + __size;
                Crispr crispr(k, _s, _e);
                crisprs.push_back(crispr);
            }
        }

        all_crisprs.insert(all_crisprs.end(), crisprs.begin(), crisprs.end());
        dyad_start += dyad_sizes[i];
    }

    done(start_time);




    return all_crisprs;
}



//        vector<unsigned int> starts(crispr_starts, crispr_starts + dyad_count);
//         printf("----begin crispr buffer-------\n");
//         for (unsigned int i = 0; i < CRISPR_BUFFER_COUNT; i++)
//         {
//             printf("%d %d/%d %d\n", k, i, CRISPR_BUFFER_COUNT, crispr_buffer[i]);
//             std::this_thread::sleep_for(std::chrono::milliseconds(1));
//         }
//         printf("------end crispr buffer-----\n");





__device__ void dyad_discovery_single_index(const char* genome, size_t genome_len, unsigned int d_index, unsigned int* dyad_buffer)
{
    for (unsigned int k = K_START; k < K_END; k++)
    {
        if (d_index + k < genome_len && is_dyad(genome, d_index, k))
        {
            unsigned int k_jump = genome_len;
            unsigned int k_index = k - K_START;
            dyad_buffer[k_index * k_jump + d_index] = d_index;
        }
    }
}


__global__ void dyad_discovery(const char* genome, size_t genome_len, unsigned int* dyad_buffer)
{
    unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;
    for (unsigned int d_index = thread_id; d_index < genome_len; d_index += stride)
        dyad_discovery_single_index(genome, genome_len, d_index, dyad_buffer);
}



vector<unsigned int> dyad_lengths(const vector<vector<unsigned int>>& all_dyads)
{
	printf("compute dyad lengths... ");
	clock_t start = clock();
    vector<unsigned int> lengths;
	lengths.reserve(all_dyads.size());
    for (const auto& vec : all_dyads)
    {
        lengths.push_back((unsigned int) vec.size());
    }
    done(start);
	return lengths;
}

vector<vector<unsigned int>> dyad_gen(char* device_genome, size_t genome_len)
{
    size_t buffer_count = genome_len * (K_END - K_START);
    unsigned int* dyad_buffer = new unsigned int[buffer_count];
    memset(dyad_buffer, 0, buffer_count * sizeof(unsigned int));

    unsigned int* device_dyad_buffer = cpush(dyad_buffer, buffer_count);
    
    dyad_discovery KERNEL_ARGS2(16, 128) (device_genome, genome_len, device_dyad_buffer);
    cwait();

    cpull(dyad_buffer, device_dyad_buffer, buffer_count);
    cufree(device_dyad_buffer);

    printf("extract dyads... ");
    clock_t start = clock();
    vector<vector<unsigned int>> all_dyads;
    for (unsigned int k = K_START; k < K_END; k++)
    {
        unsigned int hopscotch = genome_len * (k - K_START);
        vector<unsigned int> dyads;
        for (unsigned int i = 0; i < genome_len; i++)
        {
            unsigned int hopscotch_leap = hopscotch + i;
            unsigned int dyad = dyad_buffer[hopscotch_leap];
            
            if (dyad != 0)
                dyads.push_back(dyad);
        }
        all_dyads.push_back(dyads);
    }
    done(start);

    return all_dyads;
}


void print_buffer(unsigned int total_dyad_count, unsigned int* crispr_buffer)
{    
    unsigned int count = 0;
    for (unsigned int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        if (crispr_buffer[d_index * CRISPR_BUFFER + 1] == 0)
        {
            continue;
        }
        count += 1;

        printf("%d: ", d_index);
        for (unsigned int i = 0; i < CRISPR_BUFFER; i++)
        {
            printf("%d ", crispr_buffer[d_index * CRISPR_BUFFER + i]);
        }
        printf("\n");
    }
}

vector<Crispr> prospector_main_gpu(const string& genome)
{
    clock_t start;

    char* device_genome = cpush(genome.c_str(), genome.length());

    start = clock();
    vector<vector<unsigned int>> all_dyads = dyad_gen(device_genome, genome.length());
    done(start, "dyad_gen");
    
    start = clock();
    vector<Crispr> crisprs = crispr_gen(device_genome, all_dyads);
    done(start, "crispr_gen");

    cufree(device_genome);

    return crisprs;
}


vector<Crispr> prospector_main(const string& genome)
{
    printf("genome has size %zd\n", genome.size());
    
    clock_t start;
    start = clock();
    vector<Crispr> crisprs = prospector_main_gpu(genome);
    done(start, "prospector");
    return crisprs;
}


