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


double duration(clock_t start)
{
    return (clock() - start) / (double) CLOCKS_PER_SEC;
}

void dur_format(clock_t start)
{
    cout << "completed in " << duration(start) << " seconds" << endl;
}

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
    cout << "parsing genome... ";
    clock_t start = clock();
    string genome = parse_single_seq(file_path).seq;
    dur_format(start);
    return genome;
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

__global__ void kernel(int total_dyad_count, const char* genome, size_t genome_len, int* dyads, int* buffer, int* k_map, int buffer_size)
{
    for (int d_index = blockIdx.x * blockDim.x + threadIdx.x; d_index < total_dyad_count; d_index += blockDim.x * gridDim.x)
    {
        int k = k_map[d_index];
        int dyad = dyads[d_index];
        int buffer_start = d_index * buffer_size;

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

bool dyad(const char* genome, int start, int k_size)
{
    for (int i = 0; i < DYAD_MIN; i++)
    {
        char beginning_upstream = genome[start + i];
        char end_downstream = genome[start + k_size - i - 1];
        if (beginning_upstream != complement(end_downstream))
            return false;
    }
    return true;
}

vector<int> dyad_seqs(size_t genome_len, const char* genome, int k)
{
    vector<int> seqs;
    for (int i = 0; i < genome_len; i++)
    {
        if (dyad(genome, i, k))
            seqs.push_back(i);
    }
    return seqs;
}

vector<vector<int>> dyad_seqs(size_t genome_len, const char* genome, int k_start, int k_end)
{
    cout << "computing dyads... ";
    clock_t start = clock();
    vector<vector<int>> all_seqs;
    for (int k = k_start; k < k_end; k++)
        all_seqs.push_back(dyad_seqs(genome_len, genome, k));
    dur_format(start);
    return all_seqs;
}

vector<int> dyad_lengths(vector<vector<int>> all_dyads)
{
    cout << "computing dyad lengths... ";
    clock_t start = clock();
    vector<int> lengths;
    for (auto vec : all_dyads)
        lengths.push_back((int) vec.size());
    dur_format(start);
    return lengths;
}

vector<int> flatten(vector<vector<int>> all_dyads)
{
    vector<int> flattened;
    for (auto v : all_dyads)
        flattened.insert(flattened.end(), v.begin(), v.end());
    return flattened;
}

int* crispr_buffer(int total_dyad_count, int buffer_size)
{
    size_t neg_count = total_dyad_count * buffer_size;
    int* buffer;
    cudaMallocManaged(&buffer, neg_count * sizeof(int));
    for (int i = 0; i < neg_count; i++)
        buffer[i] = -1;
    return buffer;
}

bool vec_contains(vector<int> a, vector<int> b)
{
    for (auto a_elt : a)
    {
        if (find(b.begin(), b.end(), a_elt) == b.end())
            return false;
    }
    return true;
}

template <typename T>
T* push(const T* src, size_t bytes)
{
    cout << "allocating " << bytes << " bytes... ";
    clock_t start = clock();
    T* ptr = NULL;
    cudaMalloc((void**)&ptr, bytes);
    cudaMemcpy(ptr, src, bytes, cudaMemcpyHostToDevice);
    dur_format(start);
    return (T*) ptr;
}

int run(string genome_path, int min_repeats, int k_start, int k_end, int buffer_size)
{
    // HOST

    string actual_genome = parse_genome(genome_path);
    const char* genome = actual_genome.c_str();
    size_t genome_len = strlen(genome);

    vector<vector<int>> all_dyads = dyad_seqs(genome_len, genome, k_start, k_end);
    vector<int> lens = dyad_lengths(all_dyads);
    int total_dyad_count = accumulate(lens.begin(), lens.end(), 0);

    vector<int> dyads = flatten(all_dyads);
    int* buffer = crispr_buffer(total_dyad_count, buffer_size); // UNIFIED MEMORY

    vector<int> k_map;
    for (int k_index = 0; k_index < lens.size(); k_index++)
    {
        int k = k_start + k_index;
        for (int dyad_index_within_len = 0; dyad_index_within_len < lens[k_index]; dyad_index_within_len++)
        {
            k_map.push_back(k);
        }
    }

    size_t size_genome = genome_len * sizeof(char);
    size_t size_dyads = total_dyad_count * sizeof(int);
    size_t size_k_map = total_dyad_count * sizeof(int);

    // DEVICE

    char* device_genome = push(genome, size_genome);
    int* device_dyads = push(&dyads[0], size_dyads);
    int* device_k_map = push(&k_map[0], size_k_map);

    cout << "executing kernel... ";
    clock_t start = clock();
    kernel KERNEL_ARGS2(16, 1024) (total_dyad_count, device_genome, genome_len, device_dyads, buffer, device_k_map, buffer_size);
    cudaError err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
        return err;
    }
    dur_format(start);

    // HOST

    vector<vector<int>> vec_crisprs;

    // extract results
    cout << "extracting results...";
    start = clock();
    for (int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        int buffer_start = d_index * buffer_size;
        if (buffer[buffer_start + min_repeats - 1] == -1)
            continue;
        int k = k_map[d_index];

        vector<int> crispr;
        crispr.push_back(k);
        for (int i = 0; i < buffer_size; i++)
        {
            int val = buffer[buffer_start + i];
            if (val == -1)
                break;
            crispr.push_back(val);
        }
        vec_crisprs.push_back(crispr);
    }
    dur_format(start);


    cout << "pruning subset crisprs... ";
    start = clock();
    for (int i = 0; i < vec_crisprs.size(); i++)
    {
        for (int j = 0; j < vec_crisprs.size(); j++)
        {
            if (i == j)
                continue;

            if (vec_contains(vec_crisprs[i], vec_crisprs[j]))
                vec_crisprs[i][0] = -1;
        }
    }
    dur_format(start);

    cout << "results:" << endl;
    for (auto vec : vec_crisprs)
    {
        int k = vec[0];

        if (k == -1)
            continue;

        string crispr_str = "";
        crispr_str += to_string(vec[0]) + " " + to_string(vec[1]) + ":\t";
        for (auto val : vec)
            crispr_str += actual_genome.substr(val, vec[0]) + " ";

        cout << crispr_str << endl;
    }

    return 0;
}


int main()
{
    clock_t start = clock();
    
    string genome_path = R"(P:\CRISPR\data\pyogenes.fasta)";
    int min_repeats = 3;
    int k_start = 30;
    int k_end = 50;
    int buffer = 10;
    int status = run(genome_path, min_repeats, k_start, k_end, buffer);
    
    dur_format(start);

    return status;
}