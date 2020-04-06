#pragma once
#include "stdafx.h"



#define DEBUG 0



#define STOP "Z"
#define STOP_C 'Z'






const map <char, char> complement_table = {
	{'A', 'T'},
	{'T', 'A'},
	{'C', 'G'},
	{'G', 'C'},
	{'N', 'N'},
	{'n', 'n'},
};



bool subset(vector<int>, vector<int>);



template <typename T, typename Pred> vector<T> filter(const vector<T>& in, Pred predicate)
{
    auto start = time();
    vector<T> out;
    for (T el : in)
        if (predicate(el))
            out.push_back(el);
    time(start, "filter");
    return out;
}


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
    auto start = time();
    vector<T> flattened;
    for (vector<T> v : vecs) flattened.insert(flattened.end(), v.begin(), v.end());
    time(start, "flatten");
    return flattened;
}

template <typename T> bool contains(const vector<T>& target, const T& query)
{
    for (T elem : target)
        if (elem == query)
            return true;
    return false;
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

template <typename Iterable, typename Comp> void sort(Iterable& iterable, Comp comp)
{
    // auto start = time();
    sort(iterable.begin(), iterable.end(), comp);
    string str = fmt::format("sort {} items", iterable.size());
    // done(start, str);
}


map<string, string> parse_fasta(string);
string parse_fasta_single(string);

string reverse_complement(string);
int mismatch_count(string);
string seqs_to_fasta(vector<string>);
vector<string> kmerize(string, unsigned int);

namespace Util
{
    vector<string> load_genomes(string dir);
}