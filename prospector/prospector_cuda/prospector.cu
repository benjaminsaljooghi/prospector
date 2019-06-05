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
#include <numeric>

#include <optional>
#include <functional>
//
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

string parse_genome(string file_path)
{
    return parse_single_seq(file_path).seq;
}

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
    }
}

__device__ void run_analysis(int k_index, int n, const char* genome, int* buffer, int k_size)
{
    int result_index = 0;

    // Save this dyad as the beginning of a CRISPR
    buffer[k_index * BUFFER + result_index] = k_index;

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
            buffer[k_index * BUFFER + result_index] = candidate_start;

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

__global__ void kernel(int n, const char* genome, int* buffer, int k_start, int k_end)
{
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
    {
        run_analysis(i, n, genome, buffer, k_size);
    }
}


//
//
//string crispr_from_array_index(string genome, int* crisprs, int index, int k_size)
//{
//    string crispr = "";
//    for (int i = 0; i < BUFFER; i++)
//    {
//        int val = crisprs[index * BUFFER + i];
//        if (val == -1)
//        {
//            break;
//        }
//        crispr += genome.substr(val, k_size);
//        crispr += " ";
//    }
//    return crispr;
//}
//
//
//bool vec_contains(vector<int> a, vector<int> b)
//{
//    for (auto a_elt : a)
//    {
//        if (std::find(b.begin(), b.end(), a_elt) == b.end())
//        {
//            return false;
//        }
//    }
//    return true;
//}
//
//
//

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

vector<int> dyad_seqs(int genome_len, const char* genome, int k)
{
    vector<int> seqs;
    for (int i = 0; i < genome_len; i++)
    {
        if (dyad(genome, i, k))
        {
            seqs.push_back(i);
        }
    }
    return seqs;
}

vector<vector<int>> dyad_seqs(int genome_len, const char* genome, int k_start, int k_end)
{
    vector<vector<int>> all_seqs;
    for (int k = k_start; k < k_end; k++)
    {
        vector<int> seqs = dyad_seqs(genome_len, genome, k);
        all_seqs.push_back(seqs);
    }
    return all_seqs;
}

vector<int> dyad_lens(vector<vector<int>> all_dyads)
{
    vector<int> lens;
    for (auto vec : all_dyads)
    {
        lens.push_back(vec.size());
    }
    return lens;
}

int* dyads_vec_to_arr(vector<vector<int>> all_dyads)
{
    vector<int> flattened_dyads;
    for (auto vec : all_dyads)
    {
        for (auto dyad : vec)
        {
            flattened_dyads.push_back(dyad);
        }
    }
    return &flattened_dyads[0];
}

int* crispr_buffer(int total_dyad_count)
{
    size_t neg_count = total_dyad_count * BUFFER;
    int* buffer = (int*) cudaMallocManaged(neg_count * sizeof(int));
    for (int i = 0; i < neg_count; i++)
        buffer[i] = -1;
    return buffer;
}

#define K_START 35
#define K_END 40

int run()
{
    // host preprocess

    string path = R"(P:\CRISPR\data\pyogenes.fasta)";
    string actual_genome = parse_genome(path);
    const char* genome = actual_genome.c_str();
    size_t genome_len = strlen(genome);

    vector<vector<int>> all_dyads = dyad_seqs(genome_len, genome, K_START, K_END);
    vector<int> lens = dyad_lens(all_dyads);
    int total_dyad_count = accumulate(lens.begin(), lens.end(), 0);

    int* dyads = dyads_vec_to_arr(all_dyads);
    int* buffer = crispr_buffer(total_dyad_count);

    // device

    // genome
    char* device_genome = NULL;
    cudaMalloc((void**)& device_genome, genome_len);
    cudaMemcpy(device_genome, genome, genome_len, cudaMemcpyHostToDevice);

    // dyads
    int* device_dyads = NULL;
    cudaMalloc()

    // buffer
    cudaMallocManaged(&buffer, genome_len * BUFFER * sizeof(int));


    cout << "executing kernel... ";
    kernel KERNEL_ARGS2(1, 1) (genome_len, device_genome, buffer, k_size);
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
        if (buffer[k_index * BUFFER + 1] == -1)
            continue;

        // construct vec of this crispr
        vector<int> this_crispr;
        for (int i = 0; i < BUFFER; i++)
        {
            if (buffer[k_index * BUFFER + i] == -1)
            {
                break;
            }
            else
            {
                this_crispr.push_back(buffer[k_index * BUFFER + i]);
            }
        }
        vec_crisprs.push_back(this_crispr);
    }
    cout << "complete." << endl;


    //cout << "pruning subset crisprs... ";
    //for (int i = 0; i < vec_crisprs.size(); i++)
    //{
    //    for (int j = 0; j < vec_crisprs.size(); j++)
    //    {
    //        if (i == j)
    //            continue;

    //        if (vec_contains(vec_crisprs[i], vec_crisprs[j]))
    //        {
    //            vec_crisprs[i][0] = -1;
    //        }
    //    }
    //}
    //cout << "complete." << endl;


    cout << "results:" << endl;
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

int main()
{

    clock_t start;
    double duration;
    start = clock();
      
    int status = run();

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "completed in: " << duration << " seconds" << endl;

    return status;
}

