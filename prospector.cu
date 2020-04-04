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










// typedef unsigned long long ull;

typedef unsigned int ui;

#define BITS 4
#define SIZE 8

map<char, ui> scheme {
    {'A', 1},
    {'C', 2},
    {'G', 4},
    {'T', 8} 
};


ui encoded(const string& kmer)
{
    assert(kmer.size() == SIZE);
    ui e = 0;
    for (int i = 0; i < kmer.size(); i++)
        e |= scheme.at(kmer[i]) << (i * BITS);
    return e;
}

ui* encoded_genome(const string& genome)
{
    double __start = omp_get_wtime();
    ui num = genome.size() - SIZE + 1;
    ui* encoding = (ui*) malloc(sizeof(ui) * num);
    printf("genome encoding... ");
    #pragma omp parallel for
    for (ui i = 0; i < num; i++) encoding[i] = encoded(genome.substr(i, SIZE));
    done(__start);
    return encoding;
}


__host__ ui difference_cpu(const ui& _a, const ui& _b)
{
    ui _xor = _a ^ _b;
    ui _and = _xor & _a;
    ui _pc = __builtin_popcount(_and);
    return _pc;
}

__device__ unsigned char difference(const ui& _a, const ui& _b)
{
    ui _xor = _a ^ _b;
    ui _and = _xor & _a;
    unsigned char _pc = __popc(_and);
    return _pc;
}



void cwait()
{
    double start = omp_get_wtime();
    printf("kernel... ");
	cudaError err = cudaDeviceSynchronize();
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(err);
	}
    done(start);
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


__global__ void compute_query_map32(
    const ui* genome_encoding,
    const ui genome_encoding_size,
    unsigned char* query_map32
)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    // for each q, compute the mutation scores for the downrange 32 indices, packed all into a vector
    // then iterate over those scores and form the crisprs

    // do the popcounts of the initiating 8-mers have <= (8/MUTANT_TOLERANCE_RATIO=1) popcount?
    // compute the popcounts of each query compared to the next 32 indices
    // this should be an extremely efficient computation, and involves zero branching
    for (ui query = thread_id; query < genome_encoding_size - 200; query += stride)
    {
        ui q = genome_encoding[query];
        for (ui i = 0; i < 32; i++)
        {
            ui t = genome_encoding[query + K_START + SPACER_SKIP + i];
            query_map32[(query*32) + i] = difference(q, t);
        }
    }
}

// __global__ void discover_crisprs(
//         const char* genome,
//         const ui* genome_encoding,
//         const ui genome_encoding_size,
//         ui* buffer,
//         ui* starts,
//         unsigned char* sizes,
//         ui buffer_start,
//         ui query_start,
//         const ui k,
//         const ui allowed_mutations,
//         char* query_map32
//     )
// {
//     const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
//     const ui stride = blockDim.x * gridDim.x;

//     const ui thread_buffer_start = buffer_start + (thread_id * THREAD_BUFFER_COUNT);
//     const ui thread_buffer_limit = thread_buffer_start + THREAD_BUFFER_COUNT;

//     // now, we are only permitted to use those popcounts to generate queries, not allowed to look beyond the bound + 32 indices.
//     // this implementation may need to either be CPU or another kernel, slightly more complex because it involves branching

//     ui buffer_pointer = thread_buffer_start;
//     for (ui query = thread_id; query < genome_encoding_size - 200; query += stride) // -200 so we don't overflow in the inner loop
//     {
//         ui bound = query + k + SPACER_SKIP;
//         starts[query] = buffer_pointer;
//         buffer[buffer_pointer++] = query;
//         for (ui target = bound; target - bound <= SPACER_MAX; target++)
//         {
//             // for each q, compute mutation scores for the downrange 100

//             {
//                 buffer[buffer_pointer++] = target;
//                 bound = target + k + SPACER_SKIP;
//                 target = bound;
//             }
//         }

//         sizes[query] = buffer_pointer - starts[query];


//         #if DEBUG == 1
//             if (buffer_pointer  >= thread_buffer_limit)
//             {
//                 printf("FATAL: exceeded thread buffer limit. %lu/%lu\n", buffer_pointer, thread_buffer_limit);
//             }
//         #endif

//     }

// }

vector<Crispr> crispr_formation_single(const char* device_genome, const ui* device_genome_encoding, const ui genome_size, const ui genome_encoding_size, const ui k, unsigned char* device_map32)
{








    // cudaError er;

    // ui count_buffer = K_BUFFER_COUNT;
    // ui count_starts = genome_size - k + 1; // will form (at least) a singleton crispr for each i in the genome
    // ui count_sizes = genome_size - k + 1;

    // ui bytes_buffer = 4 * count_buffer;
    // ui bytes_starts = 4 * count_starts;
    // ui bytes_sizes = 1 * count_sizes;

    // // buffer info
    // unsigned int* buffer, *device_buffer;
    // unsigned int* starts, *device_starts;
    // unsigned char* sizes, *device_sizes;

    // er = cudaMallocHost(&buffer, bytes_buffer); checkCuda(er);
    // er = cudaMallocHost(&starts, bytes_starts); checkCuda(er);
    // er = cudaMallocHost(&sizes, bytes_sizes); checkCuda(er);


    // er = cudaMalloc(&device_buffer, bytes_buffer); checkCuda(er);
    // er = cudaMalloc(&device_starts, bytes_starts); checkCuda(er);
    // er = cudaMalloc(&device_sizes, bytes_sizes); checkCuda(er);


    // er = cudaMemset(device_buffer, 0, bytes_buffer); checkCuda(er);
    // er = cudaMemset(device_starts, 0, bytes_starts); checkCuda(er);
    // er = cudaMemset(device_sizes, 0, bytes_sizes); checkCuda(er);



    // // discover_crisprs KERNEL_ARGS3(C_GRID, C_BLOCK, 0) 
    // // (
    // //     device_genome,
    // //     device_genome_encoding,
    // //     genome_encoding_size,
    // //     device_buffer,
    // //     device_starts,
    // //     device_sizes,
    // //     0,
    // //     0,
    // //     k,
    // //     k / MUTANT_TOLERANCE_RATIO
    // // );

    // double start = omp_get_wtime();
    // printf("\t\tmap32... ");
    // cwait();
    // done(start);

    // // er = cudaMemcpy(buffer, device_buffer, bytes_buffer, cudaMemcpyDeviceToHost); checkCuda(er);
    // // er = cudaMemcpy(starts, device_starts, bytes_starts, cudaMemcpyDeviceToHost); checkCuda(er);
    // // er = cudaMemcpy(sizes, device_sizes, bytes_sizes, cudaMemcpyDeviceToHost); checkCuda(er);


    // start = omp_get_wtime();
    
    // printf("\t\textract...");
    // vector<Crispr> crisprs;
    // // for (ui q = 0; q < count_starts; q++)
    // // {
    // //     unsigned char __size = *(sizes + q);
    // //     if (__size >= MIN_REPEATS)
    // //     {
    // //         unsigned int start = *(starts + q);
    // //         vector<unsigned int> genome_indices;
    // //         for (int walk = 0; walk < __size; walk++)
    // //         {
    // //             genome_indices.push_back(  *(buffer + start + walk)  ); // just use indices/ptrs rather than a vector?
    // //         }
    // //         assert(__size == genome_indices.size()); // not necessary, just put it in to remind you later of double up on size in constructor
    // //         Crispr crispr(k, genome_indices, genome_indices.size()); // why the double up here..? refactor crispr?
    // //         crisprs.push_back(crispr);
    // //     }
    // // }
    // done(start);

    // return crisprs;

    // // vector<Crispr> crisprs;
    // // unsigned int __dyad_start = 0;
    // // for (unsigned int i = 0; i < K_COUNT; i++)
    // // {
    // //     unsigned int k = K_START + i;
    // //     vector<Crispr> crisprs;

    // //     for (unsigned int dyad_index = 0; dyad_index < dyad_counts[i]; dyad_index++)
    // //     {
    // //         unsigned char __size = *(sizes + __dyad_start + dyad_index);
    // //         if (__size >= MIN_REPEATS)
    // //         {
    // //             unsigned int start = *(starts + __dyad_start + dyad_index);
    // //             unsigned int* _s = buffer + start;
    // //             unsigned int* _e = _s + __size;
    // //             Crispr crispr(k, _s, _e);
    // //             crisprs.push_back(crispr);
    // //         }
    // //     }

    // //     crisprs.insert(crisprs.end(), crisprs.begin(), crisprs.end());
    // //     __dyad_start += dyad_counts[i];
    // // }


}

vector<Crispr> crispr_formation(const char* device_genome, const ui* device_genome_encoding, const ui genome_size, const ui genome_encoding_size, unsigned char* device_map32)
{
    vector<Crispr> all_crisprs;
    for (ui k = K_START; k < K_END; k++)
    {
        vector<Crispr> crisprs = crispr_formation_single(device_genome, device_genome_encoding, genome_size, genome_encoding_size, k, device_map32);
        all_crisprs.insert(all_crisprs.end(), crisprs.begin(), crisprs.end());
    }
    return all_crisprs;
}


    // for (ui query = thread_id; query < genome_encoding_size - 200; query += stride)
    // {
    //     ui q = genome_encoding[query];
    //     for (ui i = 0; i < 32; i++)
    //     {
    //         ui t = genome_encoding[query + K_START + SPACER_SKIP + i];
    //         query_map32[(query*32) + i] = difference(q, t);
    //     }
    // }

vector<ui> q_substrate(unsigned char* map32, ui genome_encoding_size)
{
    // how many qs in this map are containment oriented
    printf("q_substrate...");
    double start = omp_get_wtime();
    // ui count = 0;
    vector<ui> queries;
    for (ui query = 0; query < genome_encoding_size - 200; query++)
    {
        for (ui i = 0; i < 32; i++)
        {
            if (map32[(query*32) + i] <= 1) // 1 because 8 / 5 = 1
            {
                queries.push_back(query);
                break;
            }
        }
    }
    done(start);
    printf("%d %d\n", genome_encoding_size-200, queries.size());
    // return count;
    return queries;
}


bool mutant(const ui* genome_encoding, const ui& k, const ui& allowed_mutations, const ui& i, const ui& j)
{
    ui diff = 0;
    const ui chunks = k / SIZE;
    for (ui chunk = 0; chunk < chunks; chunk++)
    {
        ui _i = genome_encoding[i + (chunk * SIZE)];
        ui _j = genome_encoding[j + (chunk * SIZE)];
        diff += difference_cpu(_i, _j);
    }

    return diff <= (chunks * SIZE) / MUTANT_TOLERANCE_RATIO;

}


vector<Crispr> prospector_main_gpu(const string& genome)
{

    cudaDeviceReset();

    cudaError er;

    ui* genome_encoding = encoded_genome(genome);
    ui genome_encoding_size = genome.size() - SIZE + 1;

    char* device_genome;
    ui* device_genome_encoding;

    double start = omp_get_wtime();
    printf("meminit...");
        er = cudaMalloc(&device_genome, 1 * genome.length()); checkCuda(er);
        er = cudaMalloc(&device_genome_encoding, 4 * genome_encoding_size); checkCuda(er);
        er = cudaMemcpy(device_genome, genome.c_str(), 1 * genome.length(), cudaMemcpyHostToDevice); checkCuda(er);
        er = cudaMemcpy(device_genome_encoding, &genome_encoding[0], 4 * genome_encoding_size, cudaMemcpyHostToDevice); checkCuda(er);
        ui count_map32 = genome_encoding_size * 32;
        ui bytes_map32 = 1 * count_map32;
        unsigned char* map32, *device_map32;
        er = cudaMallocHost(&map32, bytes_map32); checkCuda(er);
        er = cudaMalloc(&device_map32, bytes_map32); checkCuda(er);
        er = cudaMemset(device_map32, 0, bytes_map32); checkCuda(er);
    done(start);

    compute_query_map32 KERNEL_ARGS3(C_GRID, C_BLOCK, 0)
    (
        device_genome_encoding,
        genome_encoding_size,
        device_map32
    );

    cwait();

    er = cudaMemcpy(map32, device_map32, bytes_map32, cudaMemcpyDeviceToHost); checkCuda(er);

    cudaFree(device_genome);

    vector<ui> queries = q_substrate(map32, genome_encoding_size);
    


    // for 36-mer crisprs on the cpu

    printf("crispr collection from substrate...");
    start = omp_get_wtime();
        vector<vector<ui>> crisprs;
        ui allowed_mutations = 36 / 5;
        for (ui _q = 0; _q < queries.size(); _q++)
        {
            ui q = queries[_q];

            vector<ui> crispr;
            crispr.push_back(q);

            ui bound = q + 36 + SPACER_SKIP;
            
            for (ui t = bound; t - bound <= SPACER_MAX; t++)
            {
                if (mutant(genome_encoding, 36, allowed_mutations, q, t))
                {
                    crispr.push_back(t);
                    bound = t + 36 + SPACER_SKIP;
                    t = bound;
                }
            }
            crisprs.push_back(crispr);
        }
    done(start);

    printf("crispr initialization...");
    start = omp_get_wtime();
        vector<Crispr> _crisprs;
        for (vector<ui> c : crisprs)
        {
            if (c.size() >= MIN_REPEATS)
            {
                Crispr _c(36, c, c.size());
                _crisprs.push_back(_c);   
            }
        }
    done(start);

    return _crisprs;
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


