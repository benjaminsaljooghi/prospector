#pragma once
#include "stdafx.h"


const unsigned int DEBUG_START = 0;
const unsigned int DEBUG_END = 10e7;

// const unsigned int DEBUG_START = 1085434-1;
// const unsigned int DEBUG_END = 1086855;

#define DEBUG 1


#define K_START 20
#define K_END 50
#define K_COUNT (K_END - K_START)


#define MUTANT_TOLERANCE_RATIO 5
#define REPEAT_MIN 20
#define REPEAT_MAX 60
#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10
#define REPEATS_MIN 3
#define SCAN_DOMAIN SPACER_MAX
#define ALLOW_DISCREPANT_LENGTHS false
#define MIN_REPEATS 3
//#define CRISPR_BUFFER 50
#define printf_BYTE_FORMAT_ALIGN 10


template <typename T> T mean(vector<T> scores)
{
	T sum = 0;
	for (T score : scores)
	{
		sum += score;
	}
	return sum / (double) scores.size();
}

template <typename T> vector<T> flatten(vector<vector<T>> vecs)
{
    double start = omp_get_wtime();
    vector<T> flattened;
    for (vector<T> v : vecs) flattened.insert(flattened.end(), v.begin(), v.end());
    done(start, "flatten");
    return flattened;
}

template <typename T> T most_frequent(vector<T> elements)
{
	map<T, int> frequencies;
	for (T element : elements)
		frequencies[element] += 1;

	T consensus;
	int frequency = 0;
	for (T element : elements)
	{
		if (frequencies[element] > frequency)
		{
			consensus = element;
			frequency = frequencies[element];
		}
	}
	return consensus;
}

double duration (double);
void done(double, string);
void done(double);
bool subset(vector<int>, vector<int>);