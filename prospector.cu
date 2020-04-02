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



int bits_required = 4;
int pack = 64/bits_required; // 16


map<char, unsigned long long> scheme {
    {'A', 1},
    {'C', 2},
    {'G', 4},
    {'T', 8} 
};

int ulongs_required(unsigned int k)
{
    int result = (k + pack -1) / pack;
    // printf("%d: %d\n", k, result);
    return result;
}

unsigned long long encoded(char nuc)
{
    return scheme.at(nuc);
}

typedef vector<unsigned long long> EnK;

EnK encoded(string kmer)
{
    vector<unsigned long long> encoding;
    int kmer_index = 0;
    for (int j = 0; j < ulongs_required(kmer.size()); j++)
    {
        unsigned long long e = 0;
        for (int i = 0; i < pack; i++)
        {
            e += (encoded(kmer[kmer_index++]) << (i * bits_required));
            if (kmer_index == kmer.size())
            {
                break;
            }
        }
        encoding.push_back(e);
    }
    return encoding;
}


void debug_encoding(vector<unsigned long long> encoding)
{
    for (unsigned long long l : encoding)
    {
        bitset<64> bs(l);
        cout << bs << " ";  
    }
    cout << endl;
}

vector<EnK> encoded_genome(string& genome, size_t k)
{
    double __start = omp_get_wtime();
    vector<EnK> encoding;
    for (size_t i = 0; i < genome.size() - k + 1; i++)
    {
        string kmer = genome.substr(i, k);
        encoding.push_back(encoded(kmer));
    }
    done(__start, "\t\tencoding");
    return encoding;
}

unsigned long long difference(const unsigned long long& _a, const unsigned long long& _b)
{
    unsigned long long _xor = _a ^ _b;
    unsigned long long _and = _xor & _a;
    unsigned long long _pc = __builtin_popcountll(_and);
    return _pc;
}


vector<vector<unsigned int>> substrate(unsigned int k, vector<EnK> genome_encoding)
{
    double start = omp_get_wtime();
    unsigned int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;
    int required = ulongs_required(k);
    vector<vector<unsigned int>> k_crisprs;
    for (unsigned int i = 0; i < genome_encoding.size(); i++)
    {
        vector<unsigned int> crispr; crispr.push_back(i);
        unsigned int bound = i + k + SPACER_SKIP;
        for (unsigned int j = bound; j - bound <= SPACER_MAX && j < genome_encoding.size(); j++)
        {
            int diff = difference(genome_encoding[i][0], genome_encoding[j][0]);
            if (diff > allowed_mutations)
            {
                continue;
            }

            for (int z = 1; z < required; z++)
            {
                diff += difference(genome_encoding[i][z], genome_encoding[j][z]);
            }
            if (diff <= allowed_mutations)
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

vector<Crispr> crispr_formation_pure(string genome)
{
    double start = omp_get_wtime();
    vector<vector<vector<unsigned int>>> crisprs; 
    for (unsigned int k = K_START; k < K_END; k++)
    {
        vector<EnK> genome_encoding = encoded_genome(genome, k);
        vector<vector<unsigned int>> k_crisprs = substrate(k, genome_encoding);
        crisprs.push_back(k_crisprs);

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


