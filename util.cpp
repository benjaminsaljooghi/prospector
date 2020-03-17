#include "stdafx.h"


unsigned int const DEBUG_START = 0;
unsigned int const DEBUG_END = 10e7;


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


double duration (double begin)
{
  return omp_get_wtime() - begin;
}

void done (double start, string out)
{
	printf("%s done in %.3f seconds\n", out.c_str(), duration(start));
}

void done (double start)
{
	done(start, "");
}

bool subset (vector <int> a, vector <int> b)
{
	for (int e : a)
	{
		if (find(b.begin(), b.end(), e) == b.end())
			return false;
	}
	return true;
}
