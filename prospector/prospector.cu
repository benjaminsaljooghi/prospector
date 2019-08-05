#include "stdafx.h"
#include "util.h"
#include "util.inl"

#define DYAD_MIN 5
#define REPEAT_MIN 20
#define REPEAT_MAX 60
#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10
#define REPEATS_MIN 3
#define SCAN_DOMAIN 1000
#define ALLOW_DISCREPANT_LENGTHS false
#define MIN_REPEATS 3
#define K_START 20
#define K_END 60
#define BUFFER 10

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

__global__ void discover_crisprs(int total_dyad_count, const char* genome, size_t genome_len, int* dyads, int* buffer, int* k_map, int buffer_size)
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
        if (beginning_upstream != Util::complement(end_downstream))
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


int* create_buffer(int count)
{
	int* buffer = (int*)malloc(count * sizeof(int));
	for (int i = 0; i < count; i++)
		buffer[i] = -1;
	return buffer;
}


vector<int> dyad_lengths(vector<vector<int>> all_dyads)
{
	printf("compute dyad lengths...\n");
	vector<int> lengths;
	for (auto vec : all_dyads)
		lengths.push_back((int)vec.size());
	return lengths;
}

vector<vector<int>> dyad_gen(char* device_genome, size_t genome_len, int k_start, int k_end)
{
    size_t buffer_count = genome_len * (k_end - k_start);
    int* dyad_buffer = create_buffer(buffer_count);
    
    int* device_dyad_buffer = Util::cpush(dyad_buffer, buffer_count);
    
    dyad_discovery KERNEL_ARGS2(16, 128) (device_genome, genome_len, k_start, k_end, device_dyad_buffer);
    Util::cwait();

    Util::cpull(dyad_buffer, device_dyad_buffer, buffer_count);
    Util::cfree(device_dyad_buffer);

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
    vector<int> dyads = Util::flatten(all_dyads);
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
    int* device_crispr_buffer = Util::cpush(crispr_buffer, crispr_buffer_count);
    int* device_dyads = Util::cpush(&dyads[0], total_dyad_count);
    int* device_k_map = Util::cpush(&k_map[0], total_dyad_count);

    discover_crisprs KERNEL_ARGS2(16, 128) (total_dyad_count, device_genome, genome_len, device_dyads, device_crispr_buffer, device_k_map, buffer_size);
    Util::cwait();
    
    Util::cpull(crispr_buffer, device_crispr_buffer, crispr_buffer_count);

    Util::cfree(device_crispr_buffer);
	Util::cfree(device_dyads);
	Util::cfree(device_k_map);

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

            if (Util::contains(vec_crisprs[i], vec_crisprs[j]))
                vec_crisprs[i][0] = -1;
        }
    }

    return vec_crisprs;
}

void run(string genome_path, int min_repeats, int k_start, int k_end, int buffer_size)
{
    string genome = Util::parse_fasta(genome_path).begin()->second;
    char* device_genome = Util::cpush(genome.c_str(), genome.length());

    vector<vector<int>> all_dyads = dyad_gen(device_genome, genome.length(), k_start, k_end);
    vector<vector<int>> crisps = crispr_gen(device_genome, genome.length(), k_start, k_end, min_repeats, buffer_size, all_dyads);

    Util::cfree(device_genome);

    printf("results:\n");
    for (auto vec : crisps)
    {
        int k = vec[0];
        if (k == -1)
            continue;

        string crispr_str = "";
        for (int i = 1; i < vec.size(); i++)
        {
            crispr_str += genome.substr(vec[i], k) + " ";
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
	clock_t start = clock();
	string genome_path = get_genome_path(argv);
	run(genome_path, MIN_REPEATS, K_START, K_END, BUFFER);
	printf("main completed in %.3f seconds.", Util::duration(start));
	return 0;
}