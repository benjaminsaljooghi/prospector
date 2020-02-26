#include "stdafx.h"



#define POS_STRAND true


struct Crispr
{
	int k;
	vector<int> genome_indices;
};


// general
double duration(clock_t begin);
vector<int> flatten(vector<vector<int>> vecs);
bool subset(vector<int> a, vector<int> b);


// seq
map<string, string> parse_fasta(string file_path);
char _complement(char a);
string reverse_complement(string seq);
int mismatch_count(string repeat);
string seqs_to_fasta(vector<string> seqs);
vector<string> get_kmers(string seq, int k);



// crispr
vector<string> get_repeats(Crispr crispr, string genome);
vector<string> get_spacers(Crispr crispr, string genome);
float consensus_conservation(vector<string> repeats);
float spacer_conservation(vector<string> spacers);
bool repeat_substring(Crispr crispr, int _start, int _end);
bool repeat_subset(Crispr a, Crispr b);
bool any_overlap(Crispr a, Crispr b);



// print crisprs
void print_spacers(string genome, Crispr crispr, map<string, int> spacer_scores);
void print_spacers(string genome, Crispr crispr);
void print_repeats(string genome, Crispr crispr, bool reverse_complements);
void print(string genome, vector<Crispr> crisprs);
void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores);