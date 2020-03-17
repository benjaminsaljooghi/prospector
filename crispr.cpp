

#include "stdafx.h"
#include "util.h"
#include "seq.h"
#include "crispr.h"

double bp_match_score(string a, string b, bool differential)
{
	size_t a_len = a.length();
	size_t b_len = b.length();

	size_t small_length;
	size_t big_length;

	if (a_len < b_len)
	{
		small_length = a_len;
		big_length = b_len;
	}
	else
	{
		small_length = b_len;
		big_length = a_len;
	}
	
	int matches = 0;

	for (size_t i = 0; i < small_length; i++)
		matches += a[i] == b[i] ? 1 : 0;

	int max_possible_matches = differential ? big_length : small_length;
	return (double) matches / (double) max_possible_matches;
}


double get_conservation_consensus(vector<string> repeats)
{
	string consensus = most_frequent(repeats);

	auto same_length_similarity = [] (string a, string b)
	{
		assert (a.size() == b.size());
		int matches = 0;
		for (size_t i = 0; i < a.size(); i++)
			matches += a[i] == b[i] ? 1 : 0;
		return (double) matches / (double) a.length();
	};

	vector<double> similarties;
	for (size_t i = 0; i < repeats.size(); i++)
		similarties.push_back(same_length_similarity(consensus, repeats[i]));
	
	return mean(similarties);
}


double get_conservation_spacer(vector<string> spacers)
{
	// compare each spacer against every other space but do not repeat comparisons and do not compare a spacer against itself	
	double score_sum = 0;
	int comparisons = 0;
	for (size_t i = 0; i < spacers.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			string a = spacers[i];
			string b = spacers[j];
			double score = bp_match_score(a, b, true);
			score_sum += score;
			comparisons++;
		}
	}

	double mean_score = score_sum / (double) comparisons;

	// divide the score by the number of spacers to punish having more spacers
	return mean_score / (double) spacers.size();
}


Crispr::Crispr(unsigned int _k, unsigned int* inclusive, unsigned int* exclusive)
{
	k = _k;
	genome_indices = inclusive;
	size = exclusive - inclusive;
}
		
void Crispr::update(string& genome)
{
	this->repeats = vector<string>(size);
	this->spacers = vector<string>(size-1);

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

	this->start = *genome_indices;
	this->end = (genome_indices[size-1]) + k - 1;

	this->conservation_repeats = get_conservation_consensus(repeats);
	this->conservation_spacers = get_conservation_spacer(spacers);
	this->overall_heuristic = conservation_repeats - (conservation_spacers * 2.5); // high conservation_repeats and low conservation_spacers is ideal
}

// void cache_upstream_kmers(string& genome, size_t upstream_size, unsigned int _k);

void Crispr::print(string& genome, map<string, int> spacer_scores)
{
	// header
	printf("%d - %d %d\n\n", start, end, k);
	printf("\t%fh %fr %fs\n\n", overall_heuristic, conservation_repeats, conservation_spacers);


	// repeats
	printf("\trepeats (%zd)\n", repeats.size());

	for (size_t i = 0; i < repeats.size(); i++)
	{
		string repeat = repeats[i];

		int mismatches = mismatch_count(repeat);
		int matches = repeat.length() / 2 - mismatches;
		double score = (double) matches / (double) (repeat.length() / 2);

		int start = genome_indices[i];
		int end = start + k - 1;

		int dist = i == 0 ? 0 : genome_indices[i] - (genome_indices[i-1] + k);
		
		printf("\t\t");
		printf("%d/%zd", matches, repeat.length()/2);
		printf(" %d %s %d", start, repeat.c_str(), end);
		printf(" %d", dist);
		printf(" %f\n", score);
	}
	cout << endl;
	
	
	// spacers
	printf("\tspacers (%zd)\n", spacers.size());
	
	for (string spacer : spacers)
	{
		printf("\t\t");
		printf("%d/%zd", spacer_scores[spacer], spacer.length());
		printf(" %s\n", spacer.c_str());
	}

	cout << endl;
}

bool Crispr::operator>(const Crispr& obj)
{
	return this->overall_heuristic > obj.overall_heuristic;
}

bool Crispr::operator<(const Crispr& obj)
{
	return this->overall_heuristic < obj.overall_heuristic;
}

void Crispr::cache_upstream_kmers(string genome, size_t upstream_size, unsigned int _k)
{
	string upstream = genome.substr(start - upstream_size, upstream_size);
	// this->target_kmers = get_kmers_amino(upstream, _k);

}





void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores)
{
	printf("printing %zd crisprs\n", crisprs.size());
	for (Crispr crispr : crisprs) crispr.print(genome, spacer_scores);
}


// Is the start and end of the given repeat a subset of any of the repeats of Crispr 'b'? 
bool repeat_substring(Crispr b, unsigned int start, unsigned int end)
{
	for (size_t i = 0; i < b.size; i++)
	{
		unsigned int repeat_start = b.genome_indices[i];
		unsigned int repeat_end = b.genome_indices[i] + b.k - 1;

		if (start >= repeat_start && end <= repeat_end)
		{
			return true;
		}				
	}
	return false;
}



// Are all the repeats of Crispr 'a' substrings of the repeats in Crispr 'b'?
bool repeat_subset(Crispr a, Crispr b)
{
	// a cannot be a repeat_subset of b if its k is greater than b
	if (a.k > b.k)
	{
		return false;
	}

	for (size_t i = 0; i < a.size; i++)
	{
		if (!repeat_substring(b, a.genome_indices[i], a.genome_indices[i] + a.k - 1))
		{
			return false;
		}
	}

	// all a repeats are substrings of b repeats
	return true;
}


bool any_overlap(Crispr a, Crispr b)
{

	unsigned int a_start = a.genome_indices[0];
	unsigned int a_end = a.genome_indices[a.size - 1] + a.k - 1;

	unsigned int b_start = b.genome_indices[0];
	unsigned int b_end = b.genome_indices[b.size - 1] + b.k - 1;


	bool a_before_b = a_start <= b_start;
	bool b_before_a = b_start <= a_start;

	bool a_bleeds_into_b = a_before_b && a_end >= b_start;
	bool b_bleeds_into_a = b_before_a && b_end >= a_start;

	return a_bleeds_into_b || b_bleeds_into_a;
}