#include "stdafx.h"



#define POS_STRAND true


// general
double duration(clock_t begin, clock_t end);
double duration(clock_t begin);
void done(clock_t begin, string out);
void done(clock_t begin);


vector<int> flatten(vector<vector<int>> vecs);
bool subset(vector<int> a, vector<int> b);


// seq
map<string, string> parse_fasta(string file_path);
char _complement(char a);
string reverse_complement(string seq);
int mismatch_count(string repeat);
string seqs_to_fasta(vector<string> seqs);
vector<string> get_kmers(string seq, int k);



// Crispr
float get_conservation_consensus(vector<string> repeats);
float get_conservation_spacer(vector<string> spacers);


class Crispr
{
	public:
		int k;
		vector<int> genome_indices;
		int start;
		int end;

		vector<string> repeats;
		vector<string> spacers;
		float conservation_repeats;
		float conservation_spacers;
		// float overall_heuristic; // higher the number the better

		Crispr(string genome, int _k, vector<int> _genome_indices)
		{

			k = _k;
			genome_indices = _genome_indices;

			start = genome_indices[0];
			end = genome_indices[genome_indices.size()-1] + k - 1;


			for (size_t i = 0; i < genome_indices.size(); i++)
			{
				repeats.push_back(genome.substr(genome_indices[i], k).c_str());
			}

			for (size_t i = 0; i < genome_indices.size()-1; i++)
			{
				int current_repeat_end = genome_indices[i] + k;
				int next_repeat_begin = genome_indices[i+1];
				int spacer_size = next_repeat_begin - current_repeat_end;
				spacers.push_back(genome.substr(current_repeat_end, spacer_size));
			}   
			

			conservation_repeats = get_conservation_consensus(repeats);
			conservation_spacers = get_conservation_spacer(spacers);

			// overall_heuristic = conservation_repeats - conservation_spacers; // high conservation_repeats and low conservation_spacers is ideal

		}
	private:

};

// Crispr
bool repeat_substring(Crispr crispr, int _start, int _end);
bool repeat_subset(Crispr a, Crispr b);
bool any_overlap(Crispr a, Crispr b);
bool heuristic_comparison(Crispr a, Crispr b);


// print crisprs
void print_spacers(string genome, Crispr crispr, map<string, int> spacer_scores);
void print_spacers(string genome, Crispr crispr);
void print_repeats(string genome, Crispr crispr, bool reverse_complements);
void print(string genome, vector<Crispr> crisprs);
void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores);
