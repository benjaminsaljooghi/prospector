using namespace std;


// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
//#include <thrust/device_vector.h>
#include "device_launch_parameters.h"
//#include <helper_cuda.h>

//#include "device_functions.h"

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

#include "consts.h"
#include "Sequence.h"
#include <algorithm>
#include <set>
#include <ctime>

#include <assert.h>

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
        default:
            return 'n';
    }
}

__device__ bool mutant(const char* genome, int start_a, int start_b, int k)
{
    int mutations = 0;
    int allowed_mutations = k / 10;
    for (int i = 0; i < k; i++)
    {
        if (genome[start_a + i] != genome[start_b + i] && ++mutations > allowed_mutations)
        {
            return false;
        }
    }
    return true;
}

__global__ void kernel(int total_dyad_count, const char* genome, size_t genome_len, int* dyads, int* buffer, int* k_map)
{
    for (int d_index = blockIdx.x * blockDim.x + threadIdx.x; d_index < total_dyad_count; d_index += blockDim.x * gridDim.x)
    {
        int k = k_map[d_index];
        int dyad = dyads[d_index];
        int buffer_start = d_index * BUFFER;



        int repeat_index = 0;

        // Save this dyad as the beginning of a CRISPR
        buffer[buffer_start] = dyad;

        // Search for repeats of this dyad
        int candidate = dyad + k + SPACER_SKIP;
        int countdown = SCAN_DOMAIN;
        while (countdown-- > 0 && candidate + k < genome_len)
        {
            // Is this candidate a repeat?
            if (mutant(genome, dyad, candidate, k))
            {
                // Save repeat
                repeat_index++;
                buffer[buffer_start + repeat_index] = candidate;

                // Look for the next candidate
                candidate += k + SPACER_SKIP;
                countdown = SCAN_DOMAIN;
            }
            else
            {
                candidate++;
            }
        }



    }
}

string crispr_from_buffer(string genome, int* buffer, int d_index, int k)
{
    string crispr = "";
    for (int i = 0; i < BUFFER; i++)
    {
        int val = buffer[d_index * BUFFER + i];
        if (val == -1)
        {
            break;
        }
        crispr += genome.substr(val, k);
        crispr += " ";
    }
    return crispr;
}


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

vector<int> dyad_seqs(size_t genome_len, const char* genome, int k)
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

vector<vector<int>> dyad_seqs(size_t genome_len, const char* genome, int k_start, int k_end)
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
        lens.push_back((int) vec.size());
    }
    return lens;
}

vector<int> dyads_vec_to_arr(vector<vector<int>> all_dyads)
{
    vector<int> flattened_dyads;
    for (auto vec : all_dyads)
    {
        for (auto dyad : vec)
        {
            flattened_dyads.push_back(dyad);
        }
    }
    return flattened_dyads;
}

int* crispr_buffer(int total_dyad_count)
{
    size_t neg_count = total_dyad_count * BUFFER;
    int* buffer;
    cudaMallocManaged(&buffer, neg_count * sizeof(int));
    for (int i = 0; i < neg_count; i++)
    {
        buffer[i] = -1;
    }

    return buffer;
}

#define K_START 36
#define K_END 37

int run()
{
    // HOST

    string path = R"(P:\CRISPR\data\pyogenes.fasta)";
    string actual_genome = parse_genome(path);
    const char* genome = actual_genome.c_str();
    size_t genome_len = strlen(genome);

    vector<vector<int>> all_dyads = dyad_seqs(genome_len, genome, K_START, K_END);
    vector<int> lens = dyad_lens(all_dyads);
    int total_dyad_count = accumulate(lens.begin(), lens.end(), 0);

    vector<int> dyads = dyads_vec_to_arr(all_dyads);
    int* buffer = crispr_buffer(total_dyad_count); // UNIFIED MEMORY

    vector<int> k_map;
    int dyad_index = 0;
    for (int k_index = 0; k_index < lens.size(); k_index++)
    {
        int k = K_START + k_index;
        for (int dyad_index_within_len = 0; dyad_index_within_len < lens[k_index]; dyad_index_within_len++)
        {
            k_map.push_back(k);
        }
    }

    size_t size_genome = genome_len * sizeof(char);
    size_t size_dyads = total_dyad_count * sizeof(int);
    size_t size_k_map = total_dyad_count * sizeof(int);

    // DEVICE

    // genome
    char* device_genome = NULL;
    cudaMalloc((void**)&device_genome, size_genome);
    cudaMemcpy(device_genome, genome, size_genome, cudaMemcpyHostToDevice);

    // dyads
    int* device_dyads = NULL;
    cudaMalloc((void**)&device_dyads, size_dyads);
    cudaMemcpy(device_dyads, &dyads[0], size_dyads, cudaMemcpyHostToDevice);

    // k_map
    int* device_k_map = NULL;
    cudaMalloc((void**)&device_k_map, size_k_map);
    cudaMemcpy(device_k_map, &k_map[0], size_k_map, cudaMemcpyHostToDevice);


    cout << "executing kernel... ";
    clock_t start;
    double duration;
    start = clock();

    kernel KERNEL_ARGS2(8, 512) (total_dyad_count, device_genome, genome_len, device_dyads, buffer, device_k_map);
    cudaError err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
        return err;
    }
        
    cout << "complete in " << (std::clock() - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;




    // extract results
    for (int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        int buffer_start = d_index * BUFFER;
        if (buffer[buffer_start + 1] == -1)
            continue;
        int k = k_map[d_index];
        string crispr = crispr_from_buffer(genome, buffer, d_index, k);
        cout << "crispr: " << crispr << endl << endl;
    }


    //vector<vector<int>> vec_crisprs;

    //cout << "extracting results... ";
    //for (int k_index = 0; k_index < genome_len; k_index++)
    //{
    //    if (buffer[k_index * BUFFER + 1] == -1)
    //        continue;

    //    // construct vec of this crispr
    //    vector<int> this_crispr;
    //    for (int i = 0; i < BUFFER; i++)
    //    {
    //        if (buffer[k_index * BUFFER + i] == -1)
    //        {
    //            break;
    //        }
    //        else
    //        {
    //            this_crispr.push_back(buffer[k_index * BUFFER + i]);
    //        }
    //    }
    //    vec_crisprs.push_back(this_crispr);
    //}
    //cout << "complete." << endl;


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


    //cout << "results:" << endl;
    //for (auto vec : vec_crisprs)
    //{
    //    if (vec[0] == -1)
    //        continue;

    //    cout << "crispr at: " << vec[0] << endl;

    //    for (auto val : vec)
    //    {
    //        cout << actual_genome.substr(val, k_size) << " ";
    //    }

    //    cout << endl;
    //}


    return 0;
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