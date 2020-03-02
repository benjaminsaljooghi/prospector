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
#define K_END 55
#define CRISPR_BUFFER 50
#define printf_BYTE_FORMAT_ALIGN 10


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


__device__ __host__ bool mutant(const char* genome, unsigned int start_a, unsigned int start_b, unsigned int k)
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


__device__ __host__ bool is_dyad_debug_check(unsigned int start_index)
{
    return start_index >= DEBUG_START && start_index <= DEBUG_END;
}

__device__ __host__ bool is_dyad(const char* genome, unsigned int start_index, unsigned int k)
{
    if (!is_dyad_debug_check(start_index))
    {
        return false;
    }

    unsigned int end_index = start_index + k - 1;

    unsigned int range = k/2;
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


template <typename T> void cpull(T* h, const T* d, unsigned int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;

	printf("memcpy %*zd bytes from device... ", printf_BYTE_FORMAT_ALIGN, bytes);
    clock_t start = clock();
	err = cudaMemcpy(h, d, bytes, cudaMemcpyDeviceToHost);
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
	T* ptr = NULL;

	printf("malloc %*zd bytes on device... ", printf_BYTE_FORMAT_ALIGN, bytes);
    clock_t start = clock();
	err = cudaMalloc((void**)& ptr, bytes);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to malloc device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
    done(start);

	printf("memcpy %*zd bytes to device... ", printf_BYTE_FORMAT_ALIGN, bytes);
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


// __device__ __host__ void discover_crispr(const char* genome, size_t genome_len, unsigned int k, unsigned int dyad, unsigned int* buffer, unsigned int buffer_start)
// {
//     // Search for repeats of this dyad
//     unsigned int candidate_upstream = dyad + k + SPACER_SKIP;
    
//     unsigned int repeat_index = 0;

//     // Save this dyad as the beginning of a CRISPR
//     buffer[buffer_start] = dyad;


//     unsigned int countdown = SCAN_DOMAIN;
//     while (true)
//     {

//         if (candidate_upstream + k > genome_len) 
//         {
//             break;
//         }

//         if (countdown-- == 0)
//         {
//             break;
//         }

//         // Is this candidate a repeat?
//         if (mutant(genome, dyad, candidate_upstream, k))
//         {
//             buffer[buffer_start + repeat_index++] = candidate_upstream;
//             candidate_upstream += k + SPACER_SKIP;
//             countdown = SCAN_DOMAIN;
//         }
//         else
//         {
//             candidate_upstream++;
//         }
        
//     }


//     unsigned int candidate_downstream = dyad - k - SPACER_SKIP;
//     countdown = SCAN_DOMAIN;
//     while(true)
//     {
//         if (candidate_downstream > candidate_upstream)
//         {
//             break;
//         }

//         if (countdown-- == 0)
//         {
//             break;
//         }

//         if (mutant(genome, dyad, candidate_downstream, k))
//         {
//             buffer[buffer_start + repeat_index++] = candidate_downstream;
//             candidate_downstream -= k + SPACER_SKIP;
//             countdown = SCAN_DOMAIN;
//         }
//         else
//         {
//             candidate_downstream--;
//         }
//     }

// }


// __global__ void discover_crisprs(unsigned int total_dyad_count, const char* genome, size_t genome_len, unsigned int* dyads, unsigned int* buffer, unsigned int* k_map)
// {
//     unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
//     unsigned int stride = blockDim.x * gridDim.x;

//     for (unsigned int d_index = thread_id; d_index < total_dyad_count; d_index += stride)
//     {
//         unsigned int k = k_map[d_index];
//         unsigned int dyad = dyads[d_index];
//         unsigned int buffer_start = d_index * CRISPR_BUFFER;

//         discover_crispr(genome, genome_len, k, dyad, buffer, buffer_start);
//     }
// }


__device__ void write_crispr_buffer(const char* genome, size_t genome_len, unsigned int* dyads, unsigned int query_d_index, unsigned int dyad_count, unsigned int* buffer, unsigned int buffer_start, unsigned int k)
{
    unsigned int max_distance = 100;


    unsigned int query_dyad = dyads[query_d_index];
    

    unsigned int end_of_this_crispr = query_dyad + k;

    unsigned int repeat_index = 0;

    buffer[buffer_start + repeat_index++] = query_dyad;
     
    for (int target_d_index = query_d_index + 1; target_d_index < dyad_count; target_d_index++) // this for loop goes up to the dyad count but in practice this will never happen. May want to rewrite this. while loop? or for loop up to CRISPR_BUFFER?
    {
        unsigned int target_dyad = dyads[target_d_index];
        
        // add each dyad to this array if it can be?

        
        // because I am using an unsigned int I need to check if the difference overflowed due to being WITHIN the dyad + k

        // guard against overflow
        bool within_crispr = target_dyad < end_of_this_crispr;


        unsigned int distance = target_dyad - end_of_this_crispr;

        bool too_close = distance < SPACER_SKIP; 

        // printf("%d (%d) against %d: %d\n", query_dyad, within_crispr, target_dyad, distance);
        
        if (within_crispr || too_close)
        {
            continue;
        }
        
        // if distance is less than the min AND if the query is a mutant then add.
        if (distance > max_distance)
        {
            // printf("break")
            break;
        }

        if (mutant(genome, query_dyad, target_dyad, k))
        {
            
            // printf("writing %d %d \n", query_dyad, target_dyad);
            buffer[buffer_start + repeat_index++] = target_dyad;
            end_of_this_crispr = target_dyad + k;
        }
        
    }
}

__global__ void discover_crisprs_2(const char* genome, size_t genome_len, unsigned int* dyads, unsigned int dyad_count, unsigned int* buffer, unsigned int k)
{
    unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    for (unsigned int query_d_index = thread_id; query_d_index < dyad_count; query_d_index += stride)
    {
        unsigned int buffer_start = query_d_index * CRISPR_BUFFER;
        write_crispr_buffer(genome, genome_len, dyads, query_d_index, dyad_count, buffer, buffer_start, k);
    }

}



__device__ void dyad_discovery_single_index(const char* genome, size_t genome_len, unsigned int d_index, unsigned int* dyad_buffer)
{
    for (unsigned int k = K_START; k < K_END; k++)
    {
        if (is_dyad(genome, d_index, k))
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

vector<unsigned int> dyad_lengths(vector<vector<unsigned int>> all_dyads)
{
	printf("compute dyad lengths... ");
	clock_t start = clock();
    vector<unsigned int> lengths;
	for (auto vec : all_dyads)
		lengths.push_back((unsigned int)vec.size());
    done(start);
	return lengths;
}

vector<vector<unsigned int>> dyad_gen(char* device_genome, size_t genome_len)
{
    size_t buffer_count = genome_len * (K_END - K_START);
    unsigned int* dyad_buffer = new unsigned int[buffer_count];
	fill_n(dyad_buffer, buffer_count, 0);

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



// bool comparison_routine(Crispr a, Crispr b)
// {
//     return a.genome_indices[0] < b.genome_indices[0];
// }


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



vector<Crispr> crispr_gen(string genome, char* device_genome, size_t genome_len, vector<vector<unsigned int>> all_dyads)
{
    clock_t crispr_gen_start = clock();

    clock_t start;

    vector<Crispr> all_crisprs;

    for (size_t dyad_set = 0; dyad_set < all_dyads.size(); dyad_set++)
    {     
        vector<unsigned int> dyads = all_dyads[dyad_set];
        
        printf("dyad sort..."); start = clock();
        sort(dyads.begin(), dyads.end());
        done(start);

        unsigned int k = K_START + dyad_set;

        unsigned int dyad_count = dyads.size(); printf("dyad count: %d\n", dyad_count);

        printf("buffer..."); start = clock();
        unsigned int crispr_buffer_count = dyad_count * CRISPR_BUFFER; printf("crispr buffer: %d\n", crispr_buffer_count);
        unsigned int* crispr_buffer = new unsigned int[crispr_buffer_count];
        memset(crispr_buffer, 0, crispr_buffer_count * sizeof(unsigned int));
        done(start);

        unsigned int* device_crispr_buffer = cpush(crispr_buffer, crispr_buffer_count);
        unsigned int* device_dyads = cpush(&dyads[0], dyad_count);

        discover_crisprs_2 KERNEL_ARGS2(8, 256) 
                (device_genome, genome_len, device_dyads, dyad_count, device_crispr_buffer, k);    

        cwait();
        cpull(crispr_buffer, device_crispr_buffer, crispr_buffer_count);
        cufree(device_crispr_buffer);

        vector<Crispr> k_crisprs;

        printf("extract..."); start = clock();
        for (unsigned int d_index = 0; d_index < dyad_count; d_index++)
        {

            unsigned int* __start = crispr_buffer + (d_index * CRISPR_BUFFER);

            if (*(__start + MIN_REPEATS) == 0)
                continue;

            unsigned int i;
            for (i = 0; i < CRISPR_BUFFER; i++)
            {
                if (*(__start + i) == 0) break;
            }

            // vector<unsigned int> genome_indices(__start, __start + i);
            Crispr crispr(k, __start, __start + i);
            // Crispr crispr(genome, k, genome_indices);
            k_crisprs.push_back(crispr); 
            // printf("%u %lu %lu\n", k, (unsigned long) __start, (unsigned long) __start + i);
        }

        done(start);
        printf("insert..."); start = clock();
        all_crisprs.insert(all_crisprs.end(), k_crisprs.begin(), k_crisprs.end());
        done(start);

    }

    // done(crispr_gen_start, "crispr gen");
    // exit(0);
    // return vector<Crispr>();
    return all_crisprs;
}

class Cluster
{   

    public:

        unsigned int k;
        vector<unsigned int> dyads;

        Cluster(unsigned int __k)
        {
            k = __k;
        }

        void insert(unsigned int dyad)
        {
            dyads.push_back(dyad);
            // sort(dyads.begin(), dyads.end());
        }

        unsigned int start()
        {
            return dyads[0];
        }

        unsigned int end()
        {
            return dyads[dyads.size() - 1] + k - 1;
        }
    private:

};


vector<Cluster> cluster_dyads(const char* genome, unsigned int k, vector<unsigned int> dyads)
{

    int min_distance = 100;
    sort(dyads.begin(), dyads.end());
    vector<Cluster> clusters;
    for (unsigned int dyad : dyads)
    {
        bool registered = false;
        for (size_t i = 0; i < clusters.size(); i++)
        {
            Cluster* cluster = &(clusters[i]);
            bool close_to_end = dyad > (*cluster).end() && dyad - (*cluster).end() < min_distance;
            bool close_to_start = dyad < (*cluster).start() && (*cluster).start() - dyad < min_distance; // I don't think this bool is ever true
            
            if ((close_to_end || close_to_start) && mutant(genome, (*cluster).dyads[0], dyad, k))
            {
                (*cluster).insert(dyad);
                registered = true;
                break;
            }
        }

        if (!registered)
        {
            Cluster new_cluster = Cluster(k);
            new_cluster.insert(dyad);
            clusters.push_back(new_cluster);
        }
    }

    return clusters;

}

// vector<Crispr> prospector_main_cpu(string genome)
// {
//     const char* _genome = genome.c_str();

//     unsigned int genome_len = (unsigned int) genome.size();

//     char* device_genome = cpush(genome.c_str(), genome_len);
//     vector<vector<unsigned int>> all_dyads = dyad_gen(device_genome, genome.length());

//     vector<Cluster> all_clusters;
//     for (size_t i = 0; i < all_dyads.size(); i++)
//     {
//         unsigned int k = K_START + i;

//         vector<unsigned int> dyad_set = all_dyads[i];
//         vector<Cluster> clusters = cluster_dyads(genome.c_str(), k, dyad_set);

//         for (Cluster c : clusters)
//         {
//             if (c.dyads.size() > 2)
//             {
//                 all_clusters.push_back(c);
//             }
//         }
//     }

//     vector<Crispr> crisprs;
//     for (Cluster c : all_clusters)
//     {
//         crisprs.emplace_back(genome, c.k, c.dyads);
//     }
//     return crisprs;
// }


vector<Crispr> prospector_main_gpu(string genome)
{
    // clock_t start;


    char* device_genome = cpush(genome.c_str(), genome.length());

    clock_t dyad_start = clock();
    vector<vector<unsigned int>> all_dyads = dyad_gen(device_genome, genome.length());
    clock_t dyad_end = clock();

    // vector<unsigned int> all_dyads_flat = flatten(all_dyads);

    // printf("all_dyads_flat.size(): %zd\n", all_dyads_flat.size());

    printf("dyad_gen completed in %.3f seconds.\n", duration(dyad_start, dyad_end));



    clock_t crispr_start = clock();
    vector<Crispr> crisprs = crispr_gen(genome, device_genome, genome.length(), all_dyads);
    clock_t crispr_end = clock();
    
    cufree(device_genome);

    printf("crispr_gen completed in %.3f seconds.\n", duration(crispr_start, crispr_end));

    return crisprs;

}


vector<Crispr> prospector_main(string genome)
{
    clock_t start = clock();
    vector<Crispr> crisprs = prospector_main_gpu(genome);
	printf("prospector completed in %.3f seconds.\n", duration(start));
    return crisprs;
}


