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

__device__ __host__ bool mutant(const char* genome, unsigned int start_a, unsigned int start_b, unsigned int k, unsigned int allowed_mutations)
{
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


vector<Crispr> crispr_formation_pure(const char* genome, size_t genome_len)
{
    double start = omp_get_wtime();
    vector<vector<vector<unsigned int>>> crisprs; 
    #pragma omp parallel for
    for (unsigned int k = K_START; k < K_END; k++)
    {
        unsigned int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;
        vector<vector<unsigned int>> k_crisprs;
        for (unsigned int i = 0; i < genome_len; i++)
        {
            vector<unsigned int> crispr; crispr.push_back(i);
            unsigned int bound = i + k + SPACER_SKIP;
            for (unsigned int j = bound; j - bound <= SPACER_MAX; j++)
            {
                if (mutant(genome, i, j, k, allowed_mutations))
                {
                    crispr.push_back(j);
                    bound = j + k + SPACER_SKIP;
                    j = bound;
                }
            }
            k_crisprs.push_back(crispr);
        }
        #pragma omp critical
        {
            crisprs.push_back(k_crisprs);
        }

    }
    done(start, "\t\tcrispr collection");

    start = omp_get_wtime();
    vector<Crispr> _crisprs;
    for (unsigned int i = 0; i < K_COUNT; i++)
    {
        for (vector<unsigned int> c : crisprs[i])
        {
            if (c.size() >= MIN_REPEATS)
            {
                Crispr _c(K_START + i, c, c.size());
                _crisprs.push_back(_c);   
            }
        }
    }
    done(start, "\t\tcrispr initialization");
    return _crisprs;
}


vector<Crispr> prospector_main_gpu(const string& genome)
{
    double start;
    start = omp_get_wtime();
    printf("\tcrispr formation\n");
    vector<Crispr> crisprs = crispr_formation_pure(genome.c_str(), genome.length());
    done(start, "\tcrispr formation");
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


