// CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"


#include "../util/stdafx.h"
#include "../util/util.h"
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


// #define DYAD_MIN 10
#define MISMATCH_TOLERANCE_RATIO 3
#define MUTANT_TOLERANCE_RATIO 5
#define REPEAT_MIN 20
#define REPEAT_MAX 60
#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10
#define REPEATS_MIN 3
#define SCAN_DOMAIN SPACER_MAX
#define ALLOW_DISCREPANT_LENGTHS false
#define MIN_REPEATS 3
#define K_START 20
#define K_END 60
#define BUFFER 50
#define PRINTF_BYTE_FORMAT_ALIGN 10


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


__device__ __host__ char complement(char nuc)
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




__device__ __host__ bool mutant(const char* genome, int start_a, int start_b, int k)
{
	int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;

	int mutations = 0;

	for (int i = 0; i < k; i++)
	{
        mutations += genome[start_a + i] == genome[start_b + i] ? 0 : 1;
		if (mutations > allowed_mutations)
        {
			return false;
        }
	}
	return true;
}



__device__ __host__ bool dyad(const char* genome, int start_index, int k)
{
    int mismatch_tolerance = k / MISMATCH_TOLERANCE_RATIO;

    int mismatch_count = 0;
    int end_index = start_index + k - 1;

	for (int i = 0; i < k/2; i++)
	{
		char upstream = genome[start_index + i];
		char downstream = genome[end_index - i];

        mismatch_count += upstream == complement(downstream) ? 0 : 1;

        if (mismatch_count == mismatch_tolerance)
        {
            return false;
        }
	}
	return true;
}


template <typename T> void cpull(T* h, const T* d, int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;

	printf("memcpy %*zd bytes from device... ", PRINTF_BYTE_FORMAT_ALIGN, bytes);
    clock_t start = clock();
	err = cudaMemcpy(h, d, bytes, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from device to host (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);
    
}


template <typename T> T* cpush(const T* src, int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;
	T* ptr = NULL;

	printf("malloc %*zd bytes on device... ", PRINTF_BYTE_FORMAT_ALIGN, bytes);
    clock_t start = clock();
	err = cudaMalloc((void**)& ptr, bytes);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to malloc device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);

	printf("memcpy %*zd bytes to device... ", PRINTF_BYTE_FORMAT_ALIGN, bytes);
    start = clock();
	err = cudaMemcpy(ptr, src, bytes, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);

	return (T*)ptr;
}



__global__ void discover_crisprs(int total_dyad_count, const char* genome, size_t genome_len, int* dyads, int* buffer, int* k_map)
{
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int d_index = thread_id; d_index < total_dyad_count; d_index += stride)
    {
        int k = k_map[d_index];

        // if (k != 36)
        // {
        //     continue;
        // }

        int dyad = dyads[d_index];
        int buffer_start = d_index * BUFFER;

        int repeat_index = 0;

        // Save this dyad as the beginning of a CRISPR
        buffer[buffer_start] = dyad;

        // Search for repeats of this dyad
        int candidate = dyad + k + SPACER_SKIP;


        int countdown = SCAN_DOMAIN;
        while (true)
        {

            if (candidate + k > genome_len)
            {
                break;
            }

            if (countdown-- == 0)
            {
                break;
            }


            // Is this candidate a repeat?
            if (mutant(genome, dyad, candidate, k))
            {
                // printf("saving repeat at %d\n", candidate);

                // Save repeat
                repeat_index++;
                buffer[buffer_start + repeat_index] = candidate;

                // Look for the next candidate
                candidate += k + SPACER_SKIP;
                countdown = SCAN_DOMAIN;
            }
            else
                candidate++;
        }
    }
}


__device__ void dyad_discovery_single_index(const char* genome, size_t genome_len, int d_index, int* dyad_buffer)
{
    for (int k = K_START; k < K_END; k++)
    {
        if (dyad(genome, d_index, k))
        {
            int k_jump = genome_len;
            int k_index = k - K_START;
            dyad_buffer[k_index * k_jump + d_index] = d_index;
        }
    }
}

__global__ void dyad_discovery(const char* genome, size_t genome_len, int* dyad_buffer)
{
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int d_index = thread_id; d_index < genome_len; d_index += stride)
        dyad_discovery_single_index(genome, genome_len, d_index, dyad_buffer);
}

vector<int> dyad_lengths(vector<vector<int>> all_dyads)
{
	printf("compute dyad lengths... ");
	clock_t start = clock();
    vector<int> lengths;
	for (auto vec : all_dyads)
		lengths.push_back((int)vec.size());
    done(start);
	return lengths;
}

vector<vector<int>> dyad_gen(char* device_genome, size_t genome_len)
{
    size_t buffer_count = genome_len * (K_END - K_START);
    int* dyad_buffer = new int[buffer_count];
	fill_n(dyad_buffer, buffer_count, -1);

    int* device_dyad_buffer = cpush(dyad_buffer, buffer_count);
    
    dyad_discovery KERNEL_ARGS2(16, 128) (device_genome, genome_len, device_dyad_buffer);
    cwait();

    cpull(dyad_buffer, device_dyad_buffer, buffer_count);
    cufree(device_dyad_buffer);

    printf("extract dyads... ");
    clock_t start = clock();
    vector<vector<int>> all_dyads;
    for (int k = K_START; k < K_END; k++)
    {
        int hopscotch = genome_len * (k - K_START);
        vector<int> dyads;
        for (int i = 0; i < genome_len; i++)
        {
            int hopscotch_leap = hopscotch + i;
            int dyad = dyad_buffer[hopscotch_leap];
            if (dyad != -1)
                dyads.push_back(dyad);
        }
        all_dyads.push_back(dyads);
    }
    done(start);

    return all_dyads;
}



bool comparison_routine(Crispr a, Crispr b)
{
    return a.genome_indices[0] < b.genome_indices[0];
}


void print_buffer(int total_dyad_count, int* crispr_buffer)
{    
    int count = 0;
    for (int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        if (crispr_buffer[d_index * BUFFER + 1] == -1)
        {
            continue;
        }
        count += 1;

        printf("%d: ", d_index);
        for (int i = 0; i < BUFFER; i++)
        {
            printf("%d ", crispr_buffer[d_index * BUFFER + i]);
        }
        printf("\n");
    }
}

vector<Crispr> crispr_gen(string genome, char* device_genome, size_t genome_len, vector<vector<int>> all_dyads)
{
    vector<int> lens = dyad_lengths(all_dyads);
    vector<int> dyads = flatten(all_dyads);
    int total_dyad_count = dyads.size();

    printf("generating k_map... ");
    clock_t start = clock();
    vector<int> k_map;
    for (int k_index = 0; k_index < lens.size(); k_index++)
    {
        int k = K_START + k_index;
        for (int dyad_index_within_len = 0; dyad_index_within_len < lens[k_index]; dyad_index_within_len++)
            k_map.push_back(k);
    }
    done(start);

	int crispr_buffer_count = total_dyad_count * BUFFER;
	int* crispr_buffer = new int[crispr_buffer_count];
	fill_n(crispr_buffer, crispr_buffer_count, -1);

    int* device_crispr_buffer = cpush(crispr_buffer, crispr_buffer_count);
    int* device_dyads = cpush(&dyads[0], total_dyad_count);
    int* device_k_map = cpush(&k_map[0], total_dyad_count);

    discover_crisprs KERNEL_ARGS2(8, 256) (total_dyad_count, device_genome, genome_len, device_dyads, device_crispr_buffer, device_k_map);
    cwait();
    
    cpull(crispr_buffer, device_crispr_buffer, crispr_buffer_count);

    cufree(device_crispr_buffer);
	cufree(device_dyads);
	cufree(device_k_map);



    vector<Crispr> crisprs;

    printf("extract results... ");
    start = clock();
    for (int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        // if this crispr does not have the required number of minimum repeats, then ignore it
        if (crispr_buffer[d_index * BUFFER + MIN_REPEATS] == -1)
            continue;


        vector<int> genome_indices;
        for (int i = 0; i < BUFFER; i++)
        {
            int val = crispr_buffer[d_index * BUFFER + i];
            if (val == -1)
                break;
            genome_indices.push_back(val);
        }

        Crispr crispr(genome, k_map[d_index], genome_indices);

        // Crispr.k = k_map[d_index];
        // Crispr.genome_indices = genome_indices;

        crisprs.push_back(crispr);
    
    }
	done(start);



    // printf("prune crisprs...\n");
    // vector<Crispr> pruned_crisprs;
    // for (int i = 0; i < crisprs.size(); i++)
    // {
    //     Crispr this_crispr = crisprs[i];
    //     bool this_crispr_is_a_subset = false;
    //     for (int j = 0; j < crisprs.size(); j++)
    //     {
    //         if (i == j)
    //             continue;

    //         Crispr other_crispr = crisprs[j];

    //         // if (subset(this_crispr.genome_indices, other_crispr.genome_indices))
    //         if (repeat_subset(this_crispr, other_crispr) || subset(this_crispr.genome_indices, other_crispr.genome_indices))
    //         {
    //             this_crispr_is_a_subset = true;
    //             break;
    //         }
    //     }

    //     if (!this_crispr_is_a_subset)
    //     {
    //         pruned_crisprs.push_back(this_crispr);
    //     }


    // }

    // sort(pruned_crisprs.begin(), pruned_crisprs.end(), comparison_routine);

    // return pruned_crisprs;
    
    return crisprs;
    
}


vector<Crispr> prospector_main(string genome)
{
	clock_t start = clock();

    char* device_genome = cpush(genome.c_str(), genome.length());

    clock_t dyad_start = clock();
    vector<vector<int>> all_dyads = dyad_gen(device_genome, genome.length());
    clock_t dyad_end = clock();

    printf("dyad_gen completed in %.3f seconds.\n", duration(dyad_start, dyad_end));

    clock_t crispr_start = clock();
    vector<Crispr> crisprs = crispr_gen(genome, device_genome, genome.length(), all_dyads);
    clock_t crispr_end = clock();
    
    cufree(device_genome);

    printf("crispr_gen completed in %.3f seconds.\n", duration(crispr_start, crispr_end));
	printf("prospector completed in %.3f seconds.\n", duration(start));
    return crisprs;
}


