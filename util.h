#pragma once
#include "stdafx.h"


#define DEBUG 0

#if DEBUG == 1
    #define DEBUG_START 1577778-1000
    #define DEBUG_END 1578028+1000
    // #define DEBUG_START 0
    // #define DEBUG_END 10e7
#endif






#define K_START 20
#define K_END 60
#define K_COUNT (K_END - K_START)


#define MUTANT_TOLERANCE_RATIO 6
#define REPEAT_MIN 20
#define REPEAT_MAX 60
#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP (SPACER_MIN - 1)
#define REPEATS_MIN 3
#define SCAN_DOMAIN SPACER_MAX
#define ALLOW_DISCREPANT_LENGTHS false
#define MIN_REPEATS 3
//#define CRISPR_BUFFER 50
#define printf_BYTE_FORMAT_ALIGN 10


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
void done(double, string);
void done(double);
bool subset(vector<int>, vector<int>);



template <typename T, typename Pred> vector<T> filter(vector<T> in, Pred predicate)
{
    printf("filtering %zd items...", in.size());
    double start = omp_get_wtime();
    vector<T> out;
    for (T el : in)
        if (predicate(el))
            out.push_back(el);
    done(start);
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



map<string, string> parse_fasta(string);
string parse_fasta_single(string);

string reverse_complement(string);
int mismatch_count(string);
string seqs_to_fasta(vector<string>);
vector<string> kmerize(string, unsigned int);