#include "stdafx.h"


unsigned int const DEBUG_START = 0;
unsigned int const DEBUG_END = 10e7;




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
