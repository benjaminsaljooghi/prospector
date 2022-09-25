#pragma once
#include "stdafx.h"



namespace Util
{
    using GenomeIdSequenceMap = map<string, string>;

    static const string stop = "Z";
    static const char stop_c = 'Z';

	std::vector<std::string> parse(std::string str, std::string delim);
	std::vector<std::string> parse_exact(std::string str, std::string delim);


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

	static const map<char, kmer> amino_encoding {
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


	ui difference_cpu(const kmer& _a, const kmer& _b);

	string translate_domain(const string&);
	string translate_genome(const string&, ull, ull, bool);

	bool any_overlap(ull a_start, ull a_final, ull b_start, ull b_final);

	template <typename T> void print_vector(vector<T> a)
	{
		for (ull i = 0; i < a.size(); i++)
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
		vector<T> out;
		for (T el : in)
			if (predicate(el))
				out.push_back(el);
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
		vector<T> flattened;
		for (vector<T> v : vecs) flattened.insert(flattened.end(), v.begin(), v.end());
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
		sort(iterable.begin(), iterable.end(), comp);
		string str = fmt::format("sort {} items", iterable.size());
	}

	map<string, string> load_genome(const std::filesystem::path&);
	map<string, string> parse_fasta(const string&, bool);
	//string parse_fasta_single(string);

	char complement(char nuc);
	void complement(string& seq);
	void reverse_complement(string& seq);

	
	ull mismatch_count(string);
	string seqs_to_fasta(vector<string>&);
	vector<string> kmerize(string seq, ui k);

	kmer encode_amino_kmer(const string& kmer);

	vector<kmer> encode_amino_kmers(vector<string> kmers, ull k);

    void assert_file(std::filesystem::path path);

    struct KmerBinaryRepresentation {
        uint64_t k;
        uint64_t mask;
    };

    std::string kmer_mask_init();
    KmerBinaryRepresentation kmer_to_binary(std::string &k_str);
    vector<uint64_t> kmers_to_binary(vector<string> &kmers);
}