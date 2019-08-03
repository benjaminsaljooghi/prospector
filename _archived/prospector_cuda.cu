using namespace std;


// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
//#include <thrust/device_vector.h>
#include "device_launch_parameters.h"
//#include <helper_cuda.h>

#include "device_functions.h"

#include <iostream>
#include <stdio.h>

#include <windows.h>

#include <cuda_fp16.h>
#include <fstream>
#include <string>
#include <map>

#include <vector>

#include <optional>
#include <functional>

#include "consts.h"
#include "Sequence.h"
#include <algorithm>
//#include "Crispr.h"
#include <set>
#include <ctime>
//
map<string, string> parse_fasta(string file_path)
{
    cout << "reading: " << file_path << endl;
    ifstream input(file_path);
    if (!input.good())
    {
        throw "Error opening " + file_path;
    }

    map<string, string> seqs;
    string line, name, content;
    while (getline(input, line))
    {
        if (line.empty() || line[0] == '>') // Identifier marker
        {
            if (!name.empty())
            {
                // Get what we read from the last entry
                seqs[name] = content;
                name.clear();
            }
            if (!line.empty())
            {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty())
        {
            if (line.find(' ') != string::npos) // Invalid sequence--no spaces allowed
            {
                name.clear();
                content.clear();
            }
            else
            {
                content += line;
            }
        }
    }
    if (!name.empty())
    {
        // Get what we read from the last 
        seqs[name] = content;
    }

    return seqs;
}

Sequence parse_single_seq(string file_path)
{
    map<string, string> seqs = parse_fasta(file_path);
    string seq = seqs.begin()->second;
    return Sequence(seq, 0);
}



// CUDA BEGIN

//

//
//__device__ int my_push_back(Sequence& seq)
//{
//    int insert_pt = atomicAdd(&vec_size, 1);
//    if (insert_pt >= N)
//    {
//        return -1;
//    }
//    vec[insert_pt] = seq;
//    return insert_pt;
//}

//__device__ void discover_crispr_cuda(Sequence genome, Sequence dyad)
//{
//    my_push_back(dyad);
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
//            //crispr.push_back(kmer);
//            my_push_back(kmer);
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
//            my_push_back(kmer);
//            index = kmer.start() - k - SPACER_SKIP;
//            countdown = reset;
//        }
//    }
//
//    //return crispr;
//}

constexpr int BUFFER = 20;

__device__ bool mutant(const char* genome, int start_a, int start_b, int k_size)
{ 
    int point_mutations = 0;
    int allowed = k_size / 10;
    for (int i = 0; i < k_size; i++)
    {
        if (genome[start_a + i] != genome[start_b + i] && ++point_mutations > allowed)
        {
            return false;
        }
    }
    return true;
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
    }
}



__device__ void save_repeat(int location, int repeat_index, int* crisprs)
{
    
}


__device__ void run_analysis(int k_index, int n, const char* genome, int* crisprs, int k_size)
{
    int result_index = 0;

    // Save this dyad as the beginning of a CRISPR
    crisprs[k_index * BUFFER + result_index] = k_index;

    // Search for repeats of this dyad
    int candidate_start = k_index + k_size + SPACER_SKIP;
    int countdown = SCAN_DOMAIN;
    while (countdown-- > 0)
    {
        // Guard against overflow
        if (candidate_start + k_size >= n)
            break;

        // Is this candidate a repeat?
        if (mutant(genome, k_index, candidate_start, k_size))
        {
            // Save repeat
            result_index++;
            crisprs[k_index * BUFFER + result_index] = candidate_start;

            // Look for the next candidate
            candidate_start = candidate_start + k_size + SPACER_SKIP;
            countdown = SCAN_DOMAIN;
        }
        else
        {
            candidate_start++;
        }
    }
}



__global__ void kernel(int n, const char* genome, int* crisprs, int k_size)
{
    for (int k_index = blockIdx.x * blockDim.x + threadIdx.x; k_index < n; k_index += blockDim.x * gridDim.x)
    {
        if (dyad(genome, k_index, k_size))
        {
            run_analysis(k_index, n, genome, crisprs, k_size);
        }
    }
}


string crispr_from_array_index(string genome, int* crisprs, int index, int k_size)
{
    string crispr = "";
    for (int i = 0; i < BUFFER; i++)
    {
        int val = crisprs[index * BUFFER + i];
        if (val == -1)
        {
            break;
        }
        crispr += genome.substr(val, k_size);
        crispr += " ";
    }
    return crispr;
}


bool vec_contains(vector<int> a, vector<int> b)
{
    for (auto a_elt : a)
    {
        if (std::find(b.begin(), b.end(), a_elt) == b.end())
        {
            return false;
        }
    }
    return true;
}



void for_k(int genome_len, string actual_genome, const char* genome, char* device_genome, int k_size)
{
    cout << "for k = " << k_size << endl;


    int* crisprs;
    
    cudaMallocManaged(&crisprs, genome_len * BUFFER * sizeof(int));



    for (int i = 0; i < genome_len * BUFFER; i++)
        crisprs[i] = -1;



    // kernel invoke
    int num_threads = genome_len - k_size + 1;

    cout << "executing kernel... ";
    kernel KERNEL_ARGS2(16, 1024) (genome_len, device_genome, crisprs, k_size);
    cudaError err = cudaDeviceSynchronize();
    if (cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
    }
    cout << "complete." << endl;


    vector<vector<int>> vec_crisprs;


    cout << "extracting results... ";
    for (int k_index = 0; k_index < genome_len; k_index++)
    {
        if (crisprs[k_index * BUFFER + 1] == -1)
            continue;

        // construct vec of this crispr
        vector<int> this_crispr;
        for (int i = 0; i < BUFFER; i++)
        {
            if (crisprs[k_index * BUFFER + i] == -1)
            {
                break;
            }
            else
            {
                this_crispr.push_back(crisprs[k_index * BUFFER + i]);
            }
        }
        vec_crisprs.push_back(this_crispr);
    }
    //cout << "complete." << endl;


    //cout << "pruning subset crisprs... ";
    for (int i = 0; i < vec_crisprs.size(); i++)
    {
        for (int j = 0; j < vec_crisprs.size(); j++)
        {
            if (i == j)
                continue;

            if (vec_contains(vec_crisprs[i], vec_crisprs[j]))
            {
                vec_crisprs[i][0] = -1;
            }
        }
    }
    cout << "complete." << endl;


    //cout << "results:" << endl;
    for (auto vec : vec_crisprs)
    {
        if (vec[0] == -1)
            continue;

        cout << "crispr at: " << vec[0] << endl;

        for (auto val : vec)
        {
            cout << actual_genome.substr(val, k_size) << " ";
        }

        cout << endl;
    }

    cudaFree(crisprs);

    cout << endl;
}



bool dyad(const char* genome, int start, int k_size)
{
    for (int i = 0; i < DYAD_MIN; i++)
    {
        char beginning_upstream = genome[start + i];
        char end_downstream = genome[start + k_size - i - 1];
        if (beginning_upstream != complement(end_downstream))
        {
            return false;
        }
    }
    return true;
}

vector<int> dyads(int genome_len, const char* genome, int k)
{
    vector<int> dyads;
    for (int i = 0; i < genome_len; i++)
    {
        if (dyad(genome, i, k))
        {
            dyads.push_back(i);
        }
    }
    return dyads;
}

int main()
{
    Sleep(10000);

    string path = R"(P:\CRISPR\bacteria\pyogenes.fasta)";
    Sequence seq = parse_single_seq(path);
    string actual_genome = seq.sequence();
    const char* genome = actual_genome.c_str();
    int genome_len = strlen(genome);

    clock_t start;
    double duration;
    start = clock();

    cout << endl << "allocating memory for the genome... ";
    char* device_genome = NULL;
    cudaMalloc((void**)& device_genome, genome_len);
    cudaMemcpy(device_genome, genome, genome_len,
        cudaMemcpyHostToDevice);

    cout << "complete." << endl << endl;

    for (int k = 35; k < 40; k++)
    {
        for_k(genome_len, actual_genome, genome, device_genome, k);
    }

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "printf: " << duration << '\n';
  
    return 0;
}
