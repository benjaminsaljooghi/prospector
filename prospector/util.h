// C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// C++
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include <ctime>
#include <string>
using namespace std;

// CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"

// Macros
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

#pragma region Timing

clock_t start;
#define BEGIN start = clock();
#define END printf("done in %.3f seconds\n", duration(start));

double duration(clock_t begin)
{
	return (clock() - begin) / (double)CLOCKS_PER_SEC;
}

#pragma endregion

#pragma region Fasta

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


string parse_single_seq(string file_path)
{
	map<string, string> seqs = parse_fasta(file_path);
	string seq = seqs.begin()->second;
	return seq;
}

#pragma endregion

#pragma region Vectors

vector<int> flatten(vector<vector<int>> vecs)
{
	vector<int> flattened;
	for (auto v : vecs)
		flattened.insert(flattened.end(), v.begin(), v.end());
	return flattened;
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

#pragma endregion

#pragma region CUDA


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


#pragma endregion
