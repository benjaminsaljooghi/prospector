#include "prospector.h"


#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"
#include <bitset>

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

ull difference(const ull& _a, const ull& _b)
{
    ull _xor = _a ^ _b;
    ull _and = _xor & _a;
    ull _pc = __builtin_popcountll(_and);
    return _pc;
}

bool mutant(const string& genome, const vector<ull>& genome_encoding, const unsigned int& i, const unsigned int& j, const unsigned int& k, const unsigned int& allowed_mutations)
{
    ull diff = 0;
    const int chunks = k / SIZE;
    for (int chunk = 0; chunk < chunks; chunk++)
    {
        diff += difference(genome_encoding[i + (chunk * SIZE)], genome_encoding[j + (chunk * SIZE)]);
    }

    if (diff > allowed_mutations)
    {
        return false;
    }

    
    // compute final diffs
    const int chars = k - (chunks * SIZE);
    for (int z = 0; z < chars; z++)
    {
        diff += genome[i + (chunks * SIZE) +  z] == genome[j + (chunks * SIZE) +  z] ? 0 : 1;
    }
    return diff <= allowed_mutations;
    return true;
}


// note: can repartition the 8-mer size differently depending on the current k.
// k20: 1 16-mer comparison + ...
// k21: 1 16-mer comparison + ...
// k24: 3 8-mer comparisons + 0
// k31: 3 8-mer comparisons + ...
// k32: 2 16-mer comparisons + 0
// k33: 2 16-mer comparisons + ...
// k40: 5 8-mer comparisons + 0
// k41: 4 8-mer comparisons + ...

vector<vector<unsigned int>> substrate(const string& genome, const unsigned int& k, const vector<ull>& genome_encoding)
{
    double start = omp_get_wtime();
    unsigned int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;
    vector<vector<unsigned int>> k_crisprs;
    for (unsigned int i = 0; i < genome_encoding.size(); i++)
    {
        vector<unsigned int> crispr; crispr.push_back(i);
        unsigned int bound = i + k + SPACER_SKIP;
        for (unsigned int j = bound; j - bound <= SPACER_MAX && j < genome_encoding.size(); j++)
        {
            if (mutant(genome, genome_encoding, i, j, k, allowed_mutations))
            {
                crispr.push_back(j);
                bound = j + k + SPACER_SKIP;
                j = bound;
            }
        }
        k_crisprs.push_back(crispr);
    }

    done(start, "\t\tsubstrate");
    return k_crisprs;
}

vector<Crispr> crispr_formation_pure(const string& genome)
{
    double start = omp_get_wtime();
    vector<vector<vector<unsigned int>>> crisprs; 

    vector<ull> genome_encoding = encoded_genome(genome);

    // #pragma omp parallel for
    for (unsigned int k = K_START; k < K_END; k++)
    {
        vector<vector<unsigned int>> k_crisprs = substrate(genome, k, genome_encoding);   
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
    vector<Crispr> crisprs = crispr_formation_pure(genome);
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


