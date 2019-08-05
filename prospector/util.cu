#include "stdafx.h"
#include "util.h"

double Util::duration(clock_t begin)
{
	return (clock() - begin) / (double)CLOCKS_PER_SEC;
}

map<string, string> Util::parse_fasta(string file_path)
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

vector<int> Util::flatten(vector<vector<int>> vecs)
{
	vector<int> flattened;
	for (auto v : vecs)
		flattened.insert(flattened.end(), v.begin(), v.end());
	return flattened;
}

bool Util::contains(vector<int> a, vector<int> b)
{
	for (auto a_elt : a)
	{
		if (find(b.begin(), b.end(), a_elt) == b.end())
			return false;
	}
	return true;
}

void Util::cfree(void* device_ptr)
{
	printf("executing cudafree\n");
	cudaError err = cudaFree(device_ptr);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to free device ptr (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
}

void Util::cwait()
{
	printf("waiting for kernel... ");
	clock_t start = clock();
	cudaError err = cudaDeviceSynchronize();
	printf("done in %.3f seconds\n", Util::duration(start));
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(err);
	}
}

__device__ char Util::complement(char nuc)
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


__device__ bool Util::mutant(const char* genome, int start_a, int start_b, int k)
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


__device__ bool Util::dyad(int dyad_min, const char* genome, int start, int k_size)
{
	for (int i = 0; i < dyad_min; i++)
	{
		char beginning_upstream = genome[start + i];
		char end_downstream = genome[start + k_size - i - 1];
		if (beginning_upstream != Util::complement(end_downstream))
			return false;
	}
	return true;
}