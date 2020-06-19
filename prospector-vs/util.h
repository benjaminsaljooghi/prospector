#pragma once
#include "stdafx.h"



namespace Util
{

	static const string stop = "Z";
    static const char stop_c = 'Z';

	std::vector<std::string> parse(std::string str, std::string delim);
	
	static const map <string, string> codon_table = {
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
		{"TAA", stop},
		{"TAG", stop},
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
		{"TGA", stop},
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

	static const map<char, ui> amino_encoding {
		{'F', 0},
		{'L', 1},
		{'I', 2},
		{'M', 3},
		{'V', 4},
		{'S', 5},
		{'P', 6},
		{'T', 7},
		{'A', 8},
		{'Y', 9},
		{'H', 10},
		{'Q', 11},
		{'N', 12},
		{'K', 13},
		{'D', 14},
		{'E', 15},
		{'C', 16},
		{'W', 17},
		{'R', 18},
		{'G', 19},
		//{'X', 20}, // any amino acid, not sure how to handle this
	};



	//static const map <char, char> complement_table =
	//{
	//	{'A', 'T'},
	//	{'T', 'A'},
	//	{'C', 'G'},
	//	{'G', 'C'},
	//	{'N', 'N'},
	//	{'n', 'n'},
	//};


	ui difference_cpu(const ui& _a, const ui& _b);

	string translate_domain(const string&);
	string translate_genome(const string&, ui, ui, bool);

	bool any_overlap(ui a, ui b, ui x, ui y);


	template <typename T> void print_vector(vector<T> a)
	{
		for (ui i = 0; i < a.size(); i++)
		{
			fmt::print("{}:{}\n", i, a[i]);
		}
	}

	// a is a subset of b
	template <typename T> bool subset (vector <T> a, vector <T> b)
	{
		for (int e : a)
		{
			if (find(b.begin(), b.end(), e) == b.end())
				return false;
		}
		return true;
	}

	template <typename T, typename Pred> vector<T> filter(vector<T>& in, Pred predicate)
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

	template <typename T> bool contains(const T* target_begin, ull target_size, const T& query)
	{
	    for (ull i = 0; i < target_size; i++)
		{
			if (*(target_begin+i) == query)
				return true;
		}
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

	string load_genome(string);
	//map<string, string> parse_fasta(string);
	//string parse_fasta_single(string);

	char complement(char nuc);
	void complement(string& seq);
	void reverse_complement(string& seq);

	
	int mismatch_count(string);
	string seqs_to_fasta(vector<string>);
	vector<string> kmerize(string, ui);

	ui encode_amino_kmer(const string& kmer);

	vector<ui> encode_amino_kmers(vector<string>, ui);

    //map<string, string> load_genomes(string dir);
}