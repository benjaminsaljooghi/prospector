#include "util.h"

template <typename T> void Util::cpull(T* h, const T* d, int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;

	printf("memcpy %*d bytes from device...\n", PRINTF_BYTE_FORMAT_ALIGN, (int) bytes);
	err = cudaMemcpy(h, d, bytes, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from device to host (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}
}


template <typename T> T* Util::cpush(const T* src, int count)
{
	size_t bytes = count * sizeof(T);

	cudaError err;
	T* ptr = NULL;

	printf("malloc %*d bytes on device...\n", PRINTF_BYTE_FORMAT_ALIGN, (int) bytes);
	err = cudaMalloc((void**)& ptr, bytes);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to malloc device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}

	printf("memcpy %*d bytes to device...\n", PRINTF_BYTE_FORMAT_ALIGN, (int) bytes);
	err = cudaMemcpy(ptr, src, bytes, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "failed to copy from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(err);
	}

	return (T*)ptr;
}
