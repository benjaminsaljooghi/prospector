using namespace std;


// For the CUDA runtime routines (prefixed with "cuda_")
#include "cuda.h"
#include "cuda_runtime.h"
//#include <thrust/device_vector.h>
#include "device_launch_parameters.h"
//#include <helper_cuda.h>

//#include "device_functions.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//#include <windows.h>

#include <cuda_fp16.h>
#include <fstream>
#include <string>
#include <map>

#include <vector>
#include <numeric>

#include <functional>

#include "consts.h"
#include "Sequence.h"
#include <algorithm>
#include <set>
#include <ctime>

#include <assert.h>


// Timing
clock_t start; 
#define BEGIN start = clock();
#define END printf("done in %.3f seconds\n", duration(start));

// Print formatting
#define PRINTF_BYTE_FORMAT_ALIGN 10

// Args
#define MIN_REPEATS 3
#define K_START 20
#define K_END 60
#define BUFFER 10

// CPP SAFE EXTRACT
double duration(clock_t begin)
{
    return (clock() - begin) / (double)CLOCKS_PER_SEC;
}


// CPP SAFE EXTRACT
map<string, string> parse_fasta(string file_path)
{
    cout << "reading: " << file_path << endl;
    ifstream input(file_path);
    if (!input.good())
    {
		throw runtime_error(strerror(errno));
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


// CPP SAFE EXTRACT
Sequence parse_single_seq(string file_path)
{
    map<string, string> seqs = parse_fasta(file_path);
    string seq = seqs.begin()->second;
    return Sequence(seq, 0);
}


// CPP SAFE EXTRACT
string parse_genome(string file_path)
{
    printf("parse genome...\n");
    string genome = parse_single_seq(file_path).seq;
    return genome;
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
            return false;
    }
    return true;
}

__global__ void kernel(int total_dyad_count, const char* genome, size_t genome_len, int* dyads, int* buffer, int* k_map, int buffer_size)
{
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int d_index = thread_id; d_index < total_dyad_count; d_index += stride)
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
                candidate++;
        }
    }
}

__device__ bool dyad(const char* genome, int start, int k_size)
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

__device__ void dyad_discovery_single_index(const char* genome, size_t genome_len, int d_index, int k_start, int k_end, int* dyad_buffer)
{
    for (int k = k_start; k < k_end; k++)
    {
        if (dyad(genome, d_index, k))
        {
            int k_jump = genome_len;
            int k_index = k - k_start;
            dyad_buffer[k_index * k_jump + d_index] = d_index;
        }
    }
}

__global__ void dyad_discovery(const char* genome, size_t genome_len, int k_start, int k_end, int* dyad_buffer)
{
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int d_index = thread_id; d_index < genome_len; d_index += stride)
        dyad_discovery_single_index(genome, genome_len, d_index, k_start, k_end, dyad_buffer);
}

vector<int> dyad_lengths(vector<vector<int>> all_dyads)
{
    printf("compute dyad lengths...\n");
    vector<int> lengths;
    for (auto vec : all_dyads)
        lengths.push_back((int)vec.size());
    return lengths;
}

vector<int> flatten(vector<vector<int>> all_dyads)
{
    vector<int> flattened;
    for (auto v : all_dyads)
        flattened.insert(flattened.end(), v.begin(), v.end());
    return flattened;
}


int* create_buffer(int count)
{
    int* buffer = (int*)malloc(count * sizeof(int));
    for (int i = 0; i < count; i++)
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

void cfree(void* device_ptr)
{
    printf("executing cudafree\n");
    cudaError err = cudaFree(device_ptr);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "failed to free device ptr (error code %s)!\n", cudaGetErrorString(err));
        exit(err);
    }
}

template <typename T> void pull(T* h, const T* d, int count)
{
    size_t bytes = count * sizeof(T);

    cudaError err;

    printf("memcpy %*d bytes from device...\n", PRINTF_BYTE_FORMAT_ALIGN, bytes);
    err = cudaMemcpy(h, d, bytes, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "failed to copy from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(err);
    }
}

template <typename T> T* push(const T* src, int count)
{

    size_t bytes = count * sizeof(T);

    cudaError err;
    T* ptr = NULL;
    
    printf("malloc %*d bytes on device...\n", PRINTF_BYTE_FORMAT_ALIGN, bytes);
    err = cudaMalloc((void**)& ptr, bytes);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "failed to malloc device (error code %s)!\n", cudaGetErrorString(err));
        exit(err);
    }

    printf("memcpy %*d bytes to device...\n", PRINTF_BYTE_FORMAT_ALIGN, bytes);
    err = cudaMemcpy(ptr, src, bytes, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "failed to copy from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(err);
    }

    return (T*)ptr;
}

void wait_cuda()
{
    printf("waiting for kernel... ");
    BEGIN;
    cudaError err = cudaDeviceSynchronize();
    END;
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
        exit(err);
    }
}


vector<vector<int>> dyad_gen(char* device_genome, size_t genome_len, int k_start, int k_end)
{
    size_t buffer_count = genome_len * (k_end - k_start);
    int* dyad_buffer = create_buffer(buffer_count);
    
    int* device_dyad_buffer = push(dyad_buffer, buffer_count);
    
    dyad_discovery KERNEL_ARGS2(16, 128) (device_genome, genome_len, k_start, k_end, device_dyad_buffer);
    wait_cuda();

    pull(dyad_buffer, device_dyad_buffer, buffer_count);
    cfree(device_dyad_buffer);

    printf("extract dyads...\n");
    vector<vector<int>> all_dyads;
    for (int k = k_start; k < k_end; k++)
    {
        int hopscotch = genome_len * (k - k_start);
        vector<int> dyads;
        for (int i = 0; i < genome_len; i++)
        {
            int hopscotch_leap = hopscotch + i;
            int dyad = dyad_buffer[hopscotch_leap];
            if (dyad != -1)
                dyads.push_back(dyad);
        }
        all_dyads.push_back(dyads);
    }

    return all_dyads;
}



vector<vector<int>> crispr_gen(char* device_genome, size_t genome_len, int k_start, int k_end, int min_repeats, int buffer_size, vector<vector<int>> all_dyads)
{
    vector<int> lens = dyad_lengths(all_dyads);
    vector<int> dyads = flatten(all_dyads);
    int total_dyad_count = dyads.size();

    vector<int> k_map;
    for (int k_index = 0; k_index < lens.size(); k_index++)
    {
        int k = k_start + k_index;
        for (int dyad_index_within_len = 0; dyad_index_within_len < lens[k_index]; dyad_index_within_len++)
            k_map.push_back(k);
    }


    int crispr_buffer_count = total_dyad_count * buffer_size;
    int* crispr_buffer = create_buffer(crispr_buffer_count);
    int* device_crispr_buffer = push(crispr_buffer, crispr_buffer_count);
    int* device_dyads = push(&dyads[0], total_dyad_count);
    int* device_k_map = push(&k_map[0], total_dyad_count);

    kernel KERNEL_ARGS2(16, 128) (total_dyad_count, device_genome, genome_len, device_dyads, device_crispr_buffer, device_k_map, buffer_size);
    wait_cuda();
    
    pull(crispr_buffer, device_crispr_buffer, crispr_buffer_count);

    cfree(device_crispr_buffer);
    cfree(device_dyads);
    cfree(device_k_map);

    vector<vector<int>> vec_crisprs;
    printf("extract results...\n");
    for (int d_index = 0; d_index < total_dyad_count; d_index++)
    {
        int buffer_start = d_index * buffer_size;
        if (crispr_buffer[buffer_start + min_repeats - 1] == -1)
            continue;
        int k = k_map[d_index];

        vector<int> crispr;
        crispr.push_back(k);
        for (int i = 0; i < buffer_size; i++)
        {
            int val = crispr_buffer[buffer_start + i];
            if (val == -1)
                break;
            crispr.push_back(val);
        }
        vec_crisprs.push_back(crispr);
    }
    
    printf("prune subset crisprs...\n");
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

    return vec_crisprs;
}


void run(string genome_path, int min_repeats, int k_start, int k_end, int buffer_size)
{
    string actual_genome = parse_genome(genome_path);
    const char* genome = actual_genome.c_str();
    size_t genome_len = strlen(genome);
    char* device_genome = push(genome, genome_len);

    vector<vector<int>> all_dyads = dyad_gen(device_genome, genome_len, k_start, k_end);
    vector<vector<int>> crisps = crispr_gen(device_genome, genome_len, k_start, k_end, min_repeats, buffer_size, all_dyads);

    cfree(device_genome);

    printf("results:\n");
    for (auto vec : crisps)
    {
        int k = vec[0];
        if (k == -1)
            continue;

        string crispr_str = "";
        for (int i = 1; i < vec.size(); i++)
        {
            crispr_str += actual_genome.substr(vec[i], k) + " ";
        }
        printf("%d %d:\t%s\n", k, vec[1], crispr_str.c_str());
    }

}

string get_genome_path(char** argv)
{
	string executed_program(argv[0]);
	string executed_dir = executed_program.substr(0, executed_program.find_last_of("\\/"));
	string genome_path = executed_dir + "/data/pyogenes.fasta";
	return genome_path;
}

int main(int argc, char** argv)
{
	string genome_path = get_genome_path(argv);
	clock_t start = clock();
	run(genome_path, MIN_REPEATS, K_START, K_END, BUFFER);
	printf("main completed in %.3f seconds.", duration(start));
	return 0;
}