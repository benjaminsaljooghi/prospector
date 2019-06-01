using namespace std;


// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
//#include <thrust/device_vector.h>
#include "device_launch_parameters.h"
//#include <helper_cuda.h>

#ifdef __CUDACC__
#define KERNEL_ARGS2(grid, block) <<< grid, block >>>
#define KERNEL_ARGS3(grid, block, sh_mem) <<< grid, block, sh_mem >>>
#define KERNEL_ARGS4(grid, block, sh_mem, stream) <<< grid, block, sh_mem, stream >>>
#else
#define KERNEL_ARGS2(grid, block)
#define KERNEL_ARGS3(grid, block, sh_mem)
#define KERNEL_ARGS4(grid, block, sh_mem, stream)
#endif
//
//
//
#include <iostream>
#include <stdio.h>
//#include <fstream>
//#include <string>
//#include <map>
//#include <vector>
//#include <algorithm>
//#include <optional>
//#include <functional>
//
//#include "consts.h"
//#include "Sequence.h"
//#include "Crispr.h"
//#include <set>
//
//map<string, string> parse_fasta(string file_path)
//{
//    cout << "reading: " << file_path << endl;
//    ifstream input(file_path);
//    if (!input.good())
//    {
//        throw "Error opening " + file_path;
//    }
//
//    map<string, string> seqs;
//    string line, name, content;
//    while (getline(input, line))
//    {
//        if (line.empty() || line[0] == '>') // Identifier marker
//        {
//            if (!name.empty())
//            {
//                // Get what we read from the last entry
//                seqs[name] = content;
//                name.clear();
//            }
//            if (!line.empty())
//            {
//                name = line.substr(1);
//            }
//            content.clear();
//        }
//        else if (!name.empty())
//        {
//            if (line.find(' ') != string::npos) // Invalid sequence--no spaces allowed
//            {
//                name.clear();
//                content.clear();
//            }
//            else
//            {
//                content += line;
//            }
//        }
//    }
//    if (!name.empty())
//    {
//        // Get what we read from the last 
//        seqs[name] = content;
//    }
//
//    return seqs;
//}
//
//Sequence parse_single_seq(string file_path)
//{
//    map<string, string> seqs = parse_fasta(file_path);
//    string seq = seqs.begin()->second;
//    return Sequence(seq, 0);
//}
//
//vector<string> get_kmers(string sequence, int k)
//{
//    vector<string> seqs;
//    for (size_t i = 0; i < sequence.length() - k + 1; i++)
//    {
//        seqs.push_back(sequence.substr(i, k));
//    }
//    return seqs;
//}
//
//bool mutant(Sequence a, Sequence b)
//{
//    if (!ALLOW_DISCREPANT_LENGTHS && a.length() != b.length())
//    {
//        throw exception();
//    }
//
//    int len = min(a.length(), b.length());
//
//    int allowed_point_mutations = a.length() / 10;
//    int point_mutations = 0;
//
//    for (int i = 0; i < len; i++)
//    {
//        if (a[i] != b[i] && ++point_mutations > allowed_point_mutations)
//        {
//            return false;
//        }
//    }
//    return true;
//}
//
//optional<Crispr> discover_crispr(Sequence genome, Sequence dyad)
//{
//
//    Crispr crispr;
//    crispr.add_repeat(dyad);
//
//    int k = dyad.length();
//
//    // Upstream scan
//    int index = dyad.start() + k + SPACER_SKIP;
//    const int reset = SCAN_DOMAIN;
//    int countdown = reset;
//    while (countdown-- > 0)
//    {
//        if (index + k > genome.end())
//        {
//            break;
//        }
//        Sequence kmer = genome.subseq(index++, k);
//        if (mutant(dyad, kmer))
//        {
//            crispr.add_repeat(kmer);
//            index = kmer.start() + k + SPACER_SKIP;
//            countdown = reset;
//        }
//    }
//
//    // Downstream scan
//    index = dyad.start() - k - SPACER_SKIP;
//    countdown = reset;
//    while (countdown-- > 0)
//    {
//        if (index < genome.start())
//        {
//            break;
//        }
//        Sequence kmer = genome.subseq(index--, k);
//        if (mutant(dyad, kmer))
//        {
//            crispr.add_repeat(kmer);
//            index = kmer.start() - k - SPACER_SKIP;   
//            countdown = reset;
//        }
//
//    }
//
//    if (crispr.repeats.size() >= REPEATS_MIN)
//    {
//        crispr.sort_repeats();
//        return optional<Crispr>{crispr};
//    }
//
//    return nullopt;
//}
//
//set<Crispr> discover_crisprs(Sequence genome, vector<Sequence> dyads)
//{
//    set<Crispr> crisprs;
//
//    size_t num_bytes = sizeof(vector<Sequence>);
//    for (int i = 0; i < dyads.size(); i++)
//        num_bytes += sizeof(dyads[i]);
//    cout << "dyads total " << num_bytes << " bytes" << endl;
//
//    cout << "discovering CRISPRs from " << dyads.size() << " dyads." << endl;
//    for (int i = 0; i < dyads.size(); i++)
//    {
//        Sequence dyad = dyads[i];
//        cout << "\rexamining dyad " << i << "/" << dyads.size() - 1 << " with start " << dyad.start() << "/" << genome.length();
//        optional<Crispr> crispr = discover_crispr(genome, dyad);
//        if (crispr.has_value())
//        {
//            cout << " -> CRISPR discovered at consensus start " << dyad.start() << endl;
//            crisprs.insert(*crispr);
//        }
//    }
//    cout << endl;
//    return crisprs;
//}

//__device__ thrust::device_vector<Sequence> discover_crispr_cuda(Sequence genome, Sequence dyad)
//{
//    thrust::device_vector<Sequence> crispr;
//
//    crispr.push_back(dyad);
//    int k = dyad.length();
//
//    // Upstream scan
//    int index = dyad.start() + k + SPACER_SKIP;
//    const int reset = SCAN_DOMAIN;
//    int countdown = reset;
//    while (countdown-- > 0)
//    {
//        if (index + k > genome.end())
//        {
//            break;
//        }
//        Sequence kmer = genome.subseq(index++, k);
//        if (mutant(dyad, kmer))
//        {
//            crispr.push_back(kmer);
//            index = kmer.start() + k + SPACER_SKIP;
//            countdown = reset;
//        }
//    }
//
//    // Downstream scan
//    index = dyad.start() - k - SPACER_SKIP;
//    countdown = reset;
//    while (countdown-- > 0)
//    {
//        if (index < genome.start())
//        {
//            break;
//        }
//        Sequence kmer = genome.subseq(index--, k);
//        if (mutant(dyad, kmer))
//        {
//            crispr.push_back(kmer);
//            index = kmer.start() - k - SPACER_SKIP;
//            countdown = reset;
//        }
//    }
//
//    return crispr;
//}

#define N 100

__global__ void kernel()
{
    printf("block %d, thread %d\n", blockIdx.x, threadIdx.x);
}


int main()
{
    kernel KERNEL_ARGS2(2, 10) ();
    cudaDeviceSynchronize();
    //string test_path = R"(P:\CRISPR\test_data\test.fasta)";
    //string aureus_path = R"(P:\CRISPR\bacteria\aureus.fasta)";
    //string pyogenes_path = R"(P:\CRISPR\bacteria\pyogenes.fasta)";

    //Sequence pyogenes = parse_single_seq(pyogenes_path);

    //vector<Sequence> dyads = pyogenes.dyads(30, 40);
    //set<Crispr> crisprs = discover_crisprs(pyogenes, dyads);
   
    //cout << "discovered CRISPRs: " << endl;
    //for (auto c : crisprs)
    //{
    //    cout << c.stringification() << endl;
    //}

    //return 0;
}
