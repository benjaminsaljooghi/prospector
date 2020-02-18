#include "stdafx.h"
#include "util.h"
#include "prospector.h"


// CUDA
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


#define DYAD_MIN 5
#define REPEAT_MIN 20
#define REPEAT_MAX 60
#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10
#define REPEATS_MIN 3
#define SCAN_DOMAIN 1000
#define ALLOW_DISCREPANT_LENGTHS false
#define MIN_REPEATS 3
#define K_START 20
#define K_END 60
#define BUFFER 10
#define PRINTF_BYTE_FORMAT_ALIGN 10


void cwait()
{
	printf("waiting for kernel... ");
	clock_t start = clock();
	cudaError err = cudaDeviceSynchronize();
	printf("done in %.3f seconds\n", Util::duration(start));
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(err);
	}
}

void cufree(void* device_ptr)
{
	printf("executing cudafree\n");
	cudaError err = cudaFree(device_ptr);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to free device ptr (error code %s)!\n", cudaGetErrorString(err));
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


__device__ bool mutant(const char* genome, int start_a, int start_b, int k)
{
	int mutations = 0;
	int allowed_mutations = k / 10;
	for (int i = 0; i < k; i++)
	{
		if (genome[start_a + i] != genome[start_b + i] && ++mutations > allowed_mutations)
			return false;
	}
	return true;
}


__device__ bool dyad(int dyad_min, const char* genome, int start, int k_size)
{
	for (int i = 0; i < dyad_min; i++)
	{
		char beginning_upstream = genome[start + i];
		char end_downstream = genome[start + k_size - i - 1];
		if (beginning_upstream != complement(end_downstream))
			return false;
	}
	return true;
}


template <typename T> void cpull(T* h, const T* d, int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;

	printf("memcpy %*d bytes from device...\n", PRINTF_BYTE_FORMAT_ALIGN, (int)bytes);
	err = cudaMemcpy(h, d, bytes, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from device to host (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
}


template <typename T> T* cpush(const T* src, int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;
	T* ptr = NULL;

	printf("malloc %*d bytes on device...\n", PRINTF_BYTE_FORMAT_ALIGN, (int) bytes);
	err = cudaMalloc((void**)& ptr, bytes);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to malloc device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}

	printf("memcpy %*d bytes to device...\n", PRINTF_BYTE_FORMAT_ALIGN, (int) bytes);
	err = cudaMemcpy(ptr, src, bytes, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}

	return (T*)ptr;
}



__global__ void discover_crisprs(int total_dyad_count, const char* genome, size_t genome_len, int* dyads, int* buffer, int* k_map, int buffer_size)
{
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int d_index = thread_id; d_index < total_dyad_count; d_index += stride)
    {
        int k = k_map[d_index];
        int dyad = dyads[d_index];
        int buffer_start = d_index * buffer_size;

        int repeat_index = 0;

        // Save this dyad as the beginning of a CRISPR
        buffer[buffer_start] = dyad;

        // Search for repeats of this dyad
        int candidate = dyad + k + SPACER_SKIP;
        int countdown = SCAN_DOMAIN;
        while (countdown-- > 0 && candidate + k < genome_len)
        {
            // Is this candidate a repeat?
            if (mutant(genome, dyad, candidate, k))
            {
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


__device__ void dyad_discovery_single_index(const char* genome, size_t genome_len, int d_index, int k_start, int k_end, int* dyad_buffer)
{
    for (int k = k_start; k < k_end; k++)
    {
        if (dyad(DYAD_MIN, genome, d_index, k))
        {
            int k_jump = genome_len;
            int k_index = k - k_start;
            dyad_buffer[k_index * k_jump + d_index] = d_index;
        }
    }
}

__global__ void dyad_discovery(const char* genome, size_t genome_len, int k_start, int k_end, int* dyad_buffer)
{
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int d_index = thread_id; d_index < genome_len; d_index += stride)
        dyad_discovery_single_index(genome, genome_len, d_index, k_start, k_end, dyad_buffer);
}

vector<int> dyad_lengths(vector<vector<int>> all_dyads)
{
	printf("compute dyad lengths...\n");
	vector<int> lengths;
	for (auto vec : all_dyads)
		lengths.push_back((int)vec.size());
	return lengths;
}

vector<vector<int>> dyad_gen(char* device_genome, size_t genome_len, int k_start, int k_end)
{
    size_t buffer_count = genome_len * (k_end - k_start);
    int* dyad_buffer = new int[buffer_count];
	fill_n(dyad_buffer, buffer_count, -1);

    int* device_dyad_buffer = cpush(dyad_buffer, buffer_count);
    
    dyad_discovery KERNEL_ARGS2(16, 128) (device_genome, genome_len, k_start, k_end, device_dyad_buffer);
    cwait();

    cpull(dyad_buffer, device_dyad_buffer, buffer_count);
    cufree(device_dyad_buffer);

    printf("extract dyads...\n");
    vector<vector<int>> all_dyads;
    for (int k = k_start; k < k_end; k++)
    {
        int hopscotch = genome_len * (k - k_start);
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

    return all_dyads;
}



bool comparison_routine(Util::Locus a, Util::Locus b)
{
    return a.genome_indices[0] < b.genome_indices[0];
}



vector<Util::Locus> crispr_gen(char* device_genome, size_t genome_len, int k_start, int k_end, int min_repeats, int buffer_size, vector<vector<int>> all_dyads)
{
    vector<int> lens = dyad_lengths(all_dyads);
    vector<int> dyads = Util::flatten(all_dyads);
    int total_dyad_count = dyads.size();

    vector<int> k_map;
    for (int k_index = 0; k_index < lens.size(); k_index++)
    {
        int k = k_start + k_index;
        for (int dyad_index_within_len = 0; dyad_index_within_len < lens[k_index]; dyad_index_within_len++)
            k_map.push_back(k);
    }

	int crispr_buffer_count = total_dyad_count * buffer_size;
	int* crispr_buffer = new int[crispr_buffer_count];
	fill_n(crispr_buffer, crispr_buffer_count, -1);

    int* device_crispr_buffer = cpush(crispr_buffer, crispr_buffer_count);
    int* device_dyads = cpush(&dyads[0], total_dyad_count);
    int* device_k_map = cpush(&k_map[0], total_dyad_count);

    discover_crisprs KERNEL_ARGS2(8, 256) (total_dyad_count, device_genome, genome_len, device_dyads, device_crispr_buffer, device_k_map, buffer_size);
    cwait();
    
    cpull(crispr_buffer, device_crispr_buffer, crispr_buffer_count);

    cufree(device_crispr_buffer);
	cufree(device_dyads);
	cufree(device_k_map);

    vector<Util::Locus> loci;
    printf("extract results...\n");
    for (int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        int buffer_start = d_index * buffer_size;
        if (crispr_buffer[buffer_start + min_repeats - 1] == -1)
            continue;
        int k = k_map[d_index];

        struct Util::Locus locus;
        locus.k = k;
        for (int i = 0; i < buffer_size; i++)
        {
            int val = crispr_buffer[buffer_start + i];
            if (val == -1)
                break;
            locus.genome_indices.push_back(val);
        }
        loci.push_back(locus);
    }
    

    printf("prune crisprs...\n");
    vector<Util::Locus> pruned_loci;
    for (int i = 0; i < loci.size(); i++)
    {
        Util::Locus this_crispr = loci[i];
        bool this_crispr_is_a_subset = false;
        for (int j = 0; j < loci.size(); j++)
        {
            if (i == j)
                continue;
            Util::Locus other_crispr = loci[j];

            // if (Util::subset(this_crispr.genome_indices, other_crispr.genome_indices))
            if (Util::repeat_subset(this_crispr, other_crispr) || Util::subset(this_crispr.genome_indices, other_crispr.genome_indices))
            {
                this_crispr_is_a_subset = true;
                break;
            }
        }

        if (!this_crispr_is_a_subset)
        {
            pruned_loci.push_back(this_crispr);
        }


    }

    sort(pruned_loci.begin(), pruned_loci.end(), comparison_routine);

    return pruned_loci;
}


Util::Prospection run(string genome_path, int min_repeats, int k_start, int k_end, int buffer_size)
{
    string genome = Util::parse_fasta(genome_path).begin()->second;
    char* device_genome = cpush(genome.c_str(), genome.length());

    vector<vector<int>> all_dyads = dyad_gen(device_genome, genome.length(), k_start, k_end);
    vector<Util::Locus> crisprs = crispr_gen(device_genome, genome.length(), k_start, k_end, min_repeats, buffer_size, all_dyads);

    cufree(device_genome);

    Util::Prospection prospection;
    prospection.genome = genome;
    prospection.crisprs = crisprs; 
    return prospection;
}

Util::Prospection Prospector::prospector_main()
{
	clock_t start = clock();
    string genome_path("/home/benjamin/proj/crispr-data/pyogenes.fasta");
	Util::Prospection prospection = run(genome_path, MIN_REPEATS, K_START, K_END, BUFFER);
	printf("main completed in %.3f seconds.\n", Util::duration(start));
    return prospection;
}

// int main(int argc, char** argv)
// {
    // Prospector::prospector_main();
// }