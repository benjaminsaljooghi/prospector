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












#define BITS 2
#define SIZE 16

map<char, ui> scheme {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3} 
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


// __host__ ui difference_cpu(const ui& _a, const ui& _b)
// {
//     ui _xor = _a ^ _b;
//     ui _and = _xor & _a;
//     ui _pc = __builtin_popcount(_and);
//     return _pc;
// }

// __device__ unsigned char difference(const ui& _a, const ui& _b)
// {
//     ui _xor = _a ^ _b;
//     ui _and = _xor & _a;
//     unsigned char _pc = __popc(_and);
//     return _pc;
// }


__host__ ui difference_cpu(const ui& _a, const ui& _b)
{
    ui _xor = (_a ^ _b);
    ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ui oddBits = _xor & 0x5555555555555555ull;
    ui comp = (evenBits >> 1) | oddBits;
    return __builtin_popcount(comp);
}

__device__ unsigned char difference(const ui& _a, const ui& _b)
{
    ui _xor = (_a ^ _b);
    ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ui oddBits = _xor & 0x5555555555555555ull;
    ui comp = (evenBits >> 1) | oddBits;
    return __popc(comp);
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


__global__ void compute_query_map64(
    const ui* genome_encoding,
    const ui genome_encoding_size,
    unsigned char* query_map64
)
{
    const ui thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const ui stride = blockDim.x * gridDim.x;

    // for each q, compute the mutation scores for the downrange 64 indices, packed all into a vector
    // then iterate over those scores and form the crisprs

    // do the popcounts of the initiating 8-mers have <= (8/MUTANT_TOLERANCE_RATIO=1) popcount?
    // compute the popcounts of each query compared to the next 64 indices
    // this should be an extremely efficient computation, and involves zero branching
    for (ui query = thread_id; query < genome_encoding_size - 200; query += stride)
    {
        ui q = genome_encoding[query];
        for (ui i = 0; i < 64; i++)
        {
            ui t = genome_encoding[query + K_START + SPACER_SKIP + i];
            query_map64[(query*64) + i] = difference(q, t);
        }
    }
}


vector<ui> q_substrate(unsigned char* map64, ui genome_encoding_size)
{
    // how many qs in this map are containment oriented
    printf("q_substrate...");
    double start = omp_get_wtime();
    // ui count = 0;
    vector<ui> queries;
    for (ui query = 0; query < genome_encoding_size - 200; query++)
    {
        for (ui i = 0; i < 64; i++)
        {
            if (map64[(query*64) + i] <= (SIZE / MUTANT_TOLERANCE_RATIO)) // 1 because 8 / 5 = 1
            {
                queries.push_back(query);
                break;
            }
        }
    }
    done(start);
    printf("%d %zd\n", genome_encoding_size-200, queries.size());
    // return count;
    return queries;
}


bool mutant(const char* genome, const ui* genome_encoding, const ui& k, const ui& allowed_mutations, const ui& i, const ui& j)
{
    ui diff = 0;
    const ui chunks = k / SIZE;
    // may generate a lot of crisprs that are filtered later (expensive) given that SIZE is large (16) here.
    // option is to do a more accurate mutation calculation either using per-char post the chunk division
    // or to encode entire kmers up to 64 into ull's to compute the difference efficiently.
    // post k=64 we can use the 64-long ull initially, and then compute a per-char difference afterwards.

    for (ui chunk = 0; chunk < chunks; chunk++)
    {
        ui _i = genome_encoding[i + (chunk * SIZE)];
        ui _j = genome_encoding[j + (chunk * SIZE)];
        diff += difference_cpu(_i, _j);
    }

    const ui checked_so_far = (chunks * SIZE);

    for (ui __i = checked_so_far; i < k; __i++)
    {
        diff += genome[i + checked_so_far + __i] == genome[j + checked_so_far + __i] ? 0 : 1; 
    }

    return diff <= allowed_mutations;
    // return diff <= (chunks * SIZE) / MUTANT_TOLERANCE_RATIO;
}


vector<vector<ui>> single_k_from_q_substrate(const char* genome, vector<ui> queries, ui* genome_encoding, const ui& k)
{
    vector<vector<ui>> crisprs;
    ui allowed_mutations = k / MUTANT_TOLERANCE_RATIO;

    for (ui _q = 0; _q < queries.size(); _q++)
    {
        ui q = queries[_q];

        vector<ui> crispr;
        crispr.push_back(q);

        ui bound = q + k + SPACER_SKIP;
        
        for (ui t = bound; t - bound <= SPACER_MAX; t++)
        {
            if (mutant(genome, genome_encoding, k, allowed_mutations, q, t))
            {
                crispr.push_back(t);
                bound = t + k + SPACER_SKIP;
                t = bound;
            }
        }
        crisprs.push_back(crispr);
    }
    return crisprs;
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
        ui count_map64 = genome_encoding_size * 64;
        ui bytes_map64 = 1 * count_map64;
        unsigned char* map64, *device_map64;
        er = cudaMallocHost(&map64, bytes_map64); checkCuda(er);
        er = cudaMalloc(&device_map64, bytes_map64); checkCuda(er);
        er = cudaMemset(device_map64, 0, bytes_map64); checkCuda(er);
    done(start);

    compute_query_map64 KERNEL_ARGS3(C_GRID, C_BLOCK, 0)
    (
        device_genome_encoding,
        genome_encoding_size,
        device_map64
    );

    cwait();

    er = cudaMemcpy(map64, device_map64, bytes_map64, cudaMemcpyDeviceToHost); checkCuda(er);

    cudaFree(device_genome);



    // debug
    // ui query = 1283501;
    // ui q = genome_encoding[query];
    // for (ui i = 0; i < 1000; i++)
    // {
    //     ui pos = query + K_START + SPACER_SKIP + i;
    //     ui diff = difference_cpu(genome_encoding[query], genome_encoding[pos]);

    //     printf("%s %d %d\n", genome.substr(pos, SIZE).c_str(), pos, diff);
    // }




    vector<ui> queries = q_substrate(map64, genome_encoding_size);
    
    printf("a..."); start = omp_get_wtime();
        vector<Crispr> all_crisprs;
        for (ui k = K_START; k < K_END; k++)
        {
            vector<vector<ui>> crisprs = single_k_from_q_substrate(genome.c_str(), queries, genome_encoding, k);

            for (vector<ui> c : crisprs)
            {
                if (c.size() >= MIN_REPEATS)
                {
                    Crispr _c(k, c, c.size());
                    all_crisprs.push_back(_c);   
                }
            }
        }
    done(start);


    return all_crisprs;
}


vector<Crispr> Prospector::prospector_main(const string& genome)
{
    assert(K_START >= SIZE);

    printf("genome has size %zd\n", genome.size());
    
    double start = omp_get_wtime(); 
    printf("prospector\n");




    vector<Crispr> crisprs; 
    crisprs = prospector_main_gpu(genome);
    
    
    done(start, "prospector");

    return crisprs;
}


