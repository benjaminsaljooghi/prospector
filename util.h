#pragma once
#include "stdafx.h"


#define DEBUG 0



#define STOP "Z"
#define STOP_C 'Z'



const map <string, string> codon_table = {
    {"TTT", "F"},
    {"TTC", "F"},
    {"TTA", "L"},
    {"TTG", "L"},
    {"CTT", "L"},
    {"CTC", "L"},
    {"CTA", "L"},
    {"CTG", "L"},
    {"ATT", "I"},
    {"ATC", "I"},
    {"ATA", "I"},
    {"ATG", "M"},
    {"GTT", "V"},
    {"GTC", "V"},
    {"GTA", "V"},
    {"GTG", "V"},
    {"TCT", "S"},
    {"TCC", "S"},
    {"TCA", "S"},
    {"TCG", "S"},
    {"CCT", "P"},
    {"CCC", "P"},
    {"CCA", "P"},
    {"CCG", "P"},
    {"ACT", "T"},
    {"ACC", "T"},
    {"ACA", "T"},
    {"ACG", "T"},
    {"GCT", "A"},
    {"GCC", "A"},
    {"GCA", "A"},
    {"GCG", "A"},
    {"TAT", "Y"},
    {"TAC", "Y"},
    {"TAA", STOP},
    {"TAG", STOP},
    {"CAT", "H"},
    {"CAC", "H"},
    {"CAA", "Q"},
    {"CAG", "Q"},
    {"AAT", "N"},
    {"AAC", "N"},
    {"AAA", "K"},
    {"AAG", "K"},
    {"GAT", "D"},
    {"GAC", "D"},
    {"GAA", "E"},
    {"GAG", "E"},
    {"TGT", "C"},
    {"TGC", "C"},
    {"TGA", STOP},
    {"TGG", "W"},
    {"CGT", "R"},
    {"CGC", "R"},
    {"CGA", "R"},
    {"CGG", "R"},
    {"AGT", "S"},
    {"AGC", "S"},
    {"AGA", "R"},
    {"AGG", "R"},
    {"GGT", "G"},
    {"GGC", "G"},
    {"GGA", "G"},
    {"GGG", "G"}
};

const map <char, char> complement_table = {
	{'A', 'T'},
	{'T', 'A'},
	{'C', 'G'},
	{'G', 'C'},
	{'N', 'N'},
	{'n', 'n'},
};



double duration (double);
void done(double, string, string);
void done(double, string);
void done(double);
bool subset(vector<int>, vector<int>);



template <typename T, typename Pred> vector<T> filter(const vector<T>& in, Pred predicate)
{
    double start = omp_get_wtime();
    vector<T> out;
    for (T el : in)
        if (predicate(el))
            out.push_back(el);
    done(start, "filter");
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
    double start = omp_get_wtime();
    vector<T> flattened;
    for (vector<T> v : vecs) flattened.insert(flattened.end(), v.begin(), v.end());
    done(start, "flatten");
    return flattened;
}

template <typename T> bool contains(vector<T> target, T query)
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
    double start = omp_get_wtime();
    sort(iterable.begin(), iterable.end(), comp);
    string str = fmt::format("sort {} items", iterable.size());
    done(start, str);
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