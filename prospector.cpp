// #include "prospector.h"

// #include "fmt/core.h"
// #include "fmt/format.h"


// typedef unsigned long long ull;
// typedef unsigned int ui;


// #define BITS 4
// #define SIZE 16
// map<char, ull> scheme {
//     {'A', 1},
//     {'C', 2},
//     {'G', 4},
//     {'T', 8} 
// };


// ull encoded(const string& kmer)
// {
//     assert(kmer.size() == SIZE);
//     ull e = 0;
//     for (int i = 0; i < kmer.size(); i++)
//         e += scheme.at(kmer[i]) << (i * BITS);
//     return e;
// }

// ull* encoded_genome(const string& genome)
// {
//     auto __start = time();
//     ull num = genome.size() - SIZE + 1;
//     ull* encoding = (ull*) malloc(sizeof(ull) * num);
//     printf("\t\tgenome encoding... ");
//     #pragma omp parallel for
//     for (ull i = 0; i < num; i++) encoding[i] = encoded(genome.substr(i, SIZE));
//     done(__start);
//     return encoding;
// }

// ull difference(const ull& _a, const ull& _b)
// {
//     ull _xor = _a ^ _b;
//     ull _and = _xor & _a;
//     ull _pc = __builtin_popcountll(_and);
//     return _pc;
// }

// bool mutant(const string& genome, const ull* genome_encoding, const ui& i, const ui& j, const ui& k, const ui& allowed_mutations)
// {
//     ull diff = 0;
//     const int chunks = k / SIZE;
//     for (int chunk = 0; chunk < chunks; chunk++)
//     {
//         diff += difference(genome_encoding[i + (chunk * SIZE)], genome_encoding[j + (chunk * SIZE)]);
//         if (diff > allowed_mutations)
//         {
//             return false;
//         }
//     }

//     // compute final diffs
//     const int chars = k - (chunks * SIZE);
//     for (int z = 0; z < chars; z++)
//     {
//         diff += genome[i + (chunks * SIZE) +  z] == genome[j + (chunks * SIZE) +  z] ? 0 : 1;
//     }
//     return diff <= allowed_mutations;
//     return true;
// }

// vector<vector<unsigned int>> substrate(const string& genome, const unsigned int& k, const ull* genome_encoding, const ull& genome_encoding_size)
// {
//     auto start = time();
//     unsigned int allowed_mutations = k / MUTANT_TOLERANCE_RATIO;
//     vector<vector<unsigned int>> k_crisprs;
//     printf("\t\t%d...", k);
//     for (unsigned int query = 0; query < genome_encoding_size; query++)
//     {
//         vector<unsigned int> crispr; crispr.push_back(query);
//         unsigned int bound = query + k + SPACER_SKIP;
//         for (unsigned int target = bound; target - bound <= SPACER_MAX && target < genome_encoding_size; target++)
//         { 
//             if (mutant(genome, genome_encoding, query, target, k, allowed_mutations))
//             {
//                 crispr.push_back(target);
//                 bound = target + k + SPACER_SKIP;
//                 target = bound;
//             }
//         }
//         k_crisprs.push_back(crispr);
//     }

//     done(start, "substrate");
//     return k_crisprs;
// }

// vector<Crispr> crispr_formation_pure(const string& genome)
// {
//     auto start = time();
//     vector<vector<vector<unsigned int>>> crisprs; 

//     ull* genome_encoding = encoded_genome(genome);

//     for (unsigned int k = K_START; k < K_END; k++)
//     {
//         vector<vector<unsigned int>> k_crisprs = substrate(genome, k, genome_encoding, genome.size() - k + 1);   
//         crisprs.push_back(k_crisprs);
//     }
//     done(start, "\t\tcrispr collection");

//     start = time();
//     vector<Crispr> _crisprs;
//     for (unsigned int i = 0; i < K_COUNT; i++)
//     {
//         for (vector<unsigned int> c : crisprs[i])
//         {
//             if (c.size() >= MIN_REPEATS)
//             {
//                 Crispr _c(K_START + i, c, c.size());
//                 _crisprs.push_back(_c);   
//             }
//         }
//     }
//     done(start, "\t\tcrispr initialization");
//     return _crisprs;
// }


// vector<Crispr> prospector_main_gpu(const string& genome)
// {
//     auto start;
//     start = time();
//     printf("\tcrispr formation\n");
//     vector<Crispr> crisprs = crispr_formation_pure(genome);
//     done(start, "\tcrispr formation");
//     return crisprs;
// }


// vector<Crispr> Prospector::prospector_main(const string& genome)
// {
//     printf("genome has size %zd\n", genome.size());
    
//     auto start;
//     printf("prospector\n");
//     start = time();
//     vector<Crispr> crisprs = prospector_main_gpu(genome);
//     done(start, "prospector");

//     return crisprs;
// }


