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

unsigned long long encoded(char nuc)
{
    return scheme.at(nuc);
}

typedef vector<unsigned long long> EnK;

EnK encoded(string kmer)
{
    int longs_required = (kmer.size() + pack -1) / pack;

    vector<unsigned long long> encoding;
    int kmer_index = 0;
    for (int j = 0; j < longs_required; j++)
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

unsigned long long diff(EnK a, EnK b, int size)
{
    unsigned long long _diff = 0;

    for (int i = 0; i < size; i++)
    {
        unsigned long long _a = a[i];
        unsigned long long _b = b[i];
        unsigned long long _xor = _a ^ _b;
        unsigned long long _and = _xor & _a;
        // cout << _and << endl;
        _diff += __builtin_popcountll(_and);
    }
    return _diff;
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
    vector<EnK> encoding;
    for (size_t i = 0; i < genome.size() - k + 1; i++)
    {
        string kmer = genome.substr(i, k);
        encoding.push_back(encoded(kmer));
    }
    return encoding;
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

bool mutant(EnK a, EnK b, int size, unsigned int allowed_mutations)
{
    unsigned long long _diff = diff(a, b, size);
    return _diff <= allowed_mutations;
}

vector<Crispr> crispr_formation_pure(string genome)
{
    double start = omp_get_wtime();
    vector<vector<vector<unsigned int>>> crisprs; 
    // #pragma omp parallel for
    for (unsigned int k = K_START; k < K_END; k++)
    {
        // get genome encoding
        double __start = omp_get_wtime();
        vector<EnK> genome_encoding = encoded_genome(genome, k);
        done(__start, "encoding");

        unsigned int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;
        vector<vector<unsigned int>> k_crisprs;
        for (unsigned int i = 0; i < genome_encoding.size(); i++)
        {
            vector<unsigned int> crispr; crispr.push_back(i);
            unsigned int bound = i + k + SPACER_SKIP;
            for (unsigned int j = bound; j - bound <= SPACER_MAX && j < genome_encoding.size(); j++)
            {
                if (mutant(genome_encoding[i], genome_encoding[j], genome_encoding[0].size(), allowed_mutations))
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


