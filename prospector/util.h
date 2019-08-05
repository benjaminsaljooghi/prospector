#include "stdafx.h"

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

#define PRINTF_BYTE_FORMAT_ALIGN 10

namespace Util
{
	double duration(clock_t begin);

	map<string, string> parse_fasta(string file_path);

	vector<int> flatten(vector<vector<int>> vecs);

	bool contains(vector<int> a, vector<int> b);

	void cfree(void* device_ptr);

	template <typename T> void cpull(T* h, const T* d, int count);

	template <typename T> T* cpush(const T* src, int count);

	void cwait();

	__global__ char complement(char nuc);
}


