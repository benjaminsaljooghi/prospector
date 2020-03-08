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


// general
double duration(double, double);
double duration(double);
void done(double, string);
void done(double);


template <typename T> vector<T> flatten(vector<vector<T>> vecs)
{
    double start = omp_get_wtime();
	vector<T> flattened;
	for (vector<T> v : vecs) flattened.insert(flattened.end(), v.begin(), v.end());
	done(start, "flatten");
	return flattened;
}


// Crispr
double get_conservation_consensus(vector<string> repeats);
double get_conservation_spacer(vector<string> spacers);
double get_conservation_spacer_2(vector<string> spacers);


// seq
map<string, string> parse_fasta(string file_path);
char _complement(char a);
string reverse_complement(string seq);
int mismatch_count(string repeat);
string seqs_to_fasta(vector<string> seqs);
vector<string> get_kmers(string seq, unsigned int k);
vector<string> get_kmers_amino(string dna, unsigned int k);
vector<string> sixwaytranslation(string dna);



class Crispr
{
public:

    // computed by constructor
    unsigned int* genome_indices;
    size_t size;
    unsigned int k;

    // computed by update
    unsigned int start;
    unsigned int end;
    vector<string> repeats;
    vector<string> spacers;
    double conservation_repeats;
    double conservation_spacers;
    double overall_heuristic;

    // computed by update2
    vector<string> target_kmers;

    Crispr (unsigned int _k, unsigned int* inclusive, unsigned int* exclusive)
    {
        k = _k;
        genome_indices = inclusive;
        size = exclusive - inclusive;
    }

    void update(string genome)
    {
        repeats = vector<string>(size);
        spacers = vector<string>(size-1);

        for (size_t i = 0; i < size; i++)
        {
            repeats[i] = genome.substr(genome_indices[i], k).c_str();
        }

        for (size_t i = 0; i < size - 1; i++)
        {
            unsigned int current_repeat_end = genome_indices[i] + k;
            unsigned int next_repeat_begin = genome_indices[i+1];
            unsigned int spacer_size = next_repeat_begin - current_repeat_end;
            spacers[i] = genome.substr(current_repeat_end, spacer_size);
        }

        start = *genome_indices;
        end = (genome_indices[size-1]) + k - 1;

        conservation_repeats = get_conservation_consensus(repeats);
        conservation_spacers = get_conservation_spacer(spacers);
        overall_heuristic = conservation_repeats - (conservation_spacers * 2.5); // high conservation_repeats and low conservation_spacers is ideal

    }

    void update2(string genome, size_t upstream_size)
    {
        string upstream = genome.substr(start - upstream_size, upstream_size);
        target_kmers = get_kmers_amino(upstream, k);
    }


private:

};




class Profile
{
    public:
        string name;
        string seq;
        vector<string> kmers;
        Profile(string _name, string _path, int _k);
};

class ProfileExecution
{
    public:
        Profile* profile;
        Crispr* crispr;
        vector<int> ordered_positions;
        map<string, vector<int>> locations_present;
        int hits;
        int hits_possible;

        ProfileExecution(Profile*, Crispr*);
        void to_string();
};




template <typename T> T median(vector<T> scores)
{
  size_t size = scores.size();

  if (size == 0)
  {
    return 0;  // Undefined, really.
  }
  else
  {
    sort(scores.begin(), scores.end());
    if (size % 2 == 0)
    {
      return (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
    else 
    {
      return scores[size / 2];
    }
  }
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


bool subset(vector<int> a, vector<int> b);





// Crispr
bool repeat_substring(Crispr crispr, unsigned int _start, unsigned int _end);
bool repeat_subset(Crispr a, Crispr b);
bool any_overlap(Crispr a, Crispr b);
bool perfect_genome_index_subset(Crispr a, Crispr b);
bool heuristic_comparison(Crispr a, Crispr b);


// print crisprs
void print_spacers(string genome, Crispr crispr, map<string, int> spacer_scores);
void print_spacers(string genome, Crispr crispr);
void print_repeats(string genome, Crispr crispr, bool reverse_complements);
void print(string genome, vector<Crispr> crisprs);
void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores);
