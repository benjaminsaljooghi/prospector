#include "stdafx.h"
#include "util.h"


// general

double duration(clock_t begin, clock_t end)
{
	return (end - begin) / (double)CLOCKS_PER_SEC;
}

double duration(clock_t begin)
{
	return duration(begin, clock());
}

void done(clock_t begin, string out)
{
	printf("%s done in %.3f seconds\n", out.c_str(), duration(begin));
}

void done(clock_t begin)
{
	done(begin, "");
}




vector<int> flatten(vector<vector<int>> vecs)
{
	vector<int> flattened;
	for (auto v : vecs)
		flattened.insert(flattened.end(), v.begin(), v.end());
	return flattened;
}

bool subset(vector<int> a, vector<int> b)
{
	for (int e : a)
	{
		// Attempt to find e in in b. If the pointer arrives at b.end() then e is not inside b. Therefore a is not contained within b
		if (find(b.begin(), b.end(), e) == b.end())
			return false;
	}
	// In all cases, elements within a were found in b. Therefore, a is a subset of b
	return true;
}



map<char, char> complement_table;

// seq


map<string, string> parse_fasta(string file_path)
{
	cout << "reading: " << file_path << endl;
	ifstream input(file_path);
	if (!input.good())
	{
		throw runtime_error(strerror(errno));
	}

	map<string, string> seqs;
	string line, name, content;
	while (getline(input, line))
	{
		if (line.empty() || line[0] == '>') // Identifier marker
		{
			if (!name.empty())
			{
				// Get what we read from the last entry
				seqs[name] = content;
				name.clear();
			}
			if (!line.empty())
			{
				name = line.substr(1);
			}
			content.clear();
		}
		else if (!name.empty())
		{
			if (line.find(' ') != string::npos) // Invalid sequence--no spaces allowed
			{
				name.clear();
				content.clear();
			}
			if (line.find("\r") != string::npos)
			{
				throw runtime_error("File contains carriage returns (CRLF). Please reformat to LF.");
			}
			else
			{
				content += line;
			}
		}
	}
	if (!name.empty())
	{
		// Get what we read from the last 
		seqs[name] = content;
	}

	return seqs;
}




char _complement(char a)
{
	complement_table['A'] = 'T';
	complement_table['T'] = 'A';
	complement_table['C'] = 'G';
	complement_table['G'] = 'C';
	complement_table['N'] = 'N';
	complement_table['n'] = 'n';

	return complement_table[a];
}

string reverse_complement(string seq)
{

    size_t len = seq.length();
    string complement = seq;

    for (size_t i = 0; i < len; i++)
    {
        complement[i] = _complement(seq[i]);
    }

	reverse(complement.begin(), complement.end());
	return complement;
}

int mismatch_count(string repeat)
{
	int _count = 0;

	size_t k = repeat.size();
	int start_index = 0;
	int end_index = start_index + repeat.size() - 1;


	for (size_t i = 0; i < k/2; i++)
	{
		char upstream = repeat[start_index + i];
		char downstream = repeat[end_index - i];

		_count += upstream == _complement(downstream) ? 0 : 1;
	}
	return _count;
}

string seqs_to_fasta(vector<string> seqs)
{
    ostringstream string_stream;
    for (size_t i = 0; i < seqs.size(); i++)
    {
        string_stream << ">" << i << endl;
        string_stream << seqs[i] << endl;  
    }
    return string_stream.str();
}

vector<string> get_kmers(string seq, int k)
{
	vector<string> kmers;
	for (size_t i = 0; i < seq.length() - k + 1; i++)
	{
		kmers.push_back(seq.substr(i, k));
	}
	return kmers;
}



// crispr


string most_frequent(vector<string> repeats)
{
	map<string, int> frequencies;
	for (string repeat : repeats)
	{
		frequencies[repeat] += 1;
	}

	string consensus;
	int frequency = 0;
	for (string repeat : repeats)
	{
		if (frequencies[repeat] > frequency)
		{
			consensus = repeat;
			frequency = frequencies[repeat];
		}
	}

	return consensus;

}


float similarity(string a, string b)
{
	if (a.length() != b.length())
	{
		return -1;
	}
	

	int matches = 0;
	for (size_t i = 0; i < a.length(); i++)
	{
		matches += a[i] == b[i] ? 1 : 0; 
	}
	
	
	return (float) matches / (float) a.length();
}


float get_conservation_consensus(vector<string> repeats)
{
	string consensus = most_frequent(repeats);

	float similarity_sum = 0;
	for (size_t i = 0; i < repeats.size(); i++)
	{
		similarity_sum += similarity(consensus, repeats[i]);

	}

	return similarity_sum / (float) repeats.size();
}

float differential_length_string_comparison(string a, string b)
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
	{
		matches += a[i] == b[i] ? 1 : 0;
	}

	return (float) matches / (float) big_length;
	
}

float get_conservation_spacer(vector<string> spacers)
{

	// float similarity_mean_sum = 0;

	// for (int i = 0; i < spacers.size(); i++)
	// {
	// 	// what is the mean similarity of this spacer against all other spacers?
	// 	float similarity_sum = 0;
	// 	for (int j = 0; j < spacers.size(); j++)
	// 	{
	// 		if (i == j)
	// 		{
	// 			continue;
	// 		}
	// 		string a = spacers[i];
	// 		string b = spacers[j];

	// 		similarity_sum += differential_length_similarity(a, b);
	// 	}

	// 	// what was the mean similarity for this spacer against all other spacers?
	// 	float similarity_mean = similarity_sum / ((float) spacers.size() - 1);

	// 	// add that to the mean
	// 	similarity_mean_sum += similarity_mean;
	// }		

	// // what was the overall mean similarity?
	// return similarity_mean_sum / (float) spacers.size(); 




	// we want this to be calculcated such that if there are more spacers then the entropy increases, reducing the conservation score (which is good)
	


	// compare each spacer against every other space but do not repeat comparisons and do not compare a spacer against itself
	
	float score_sum = 0;
	int comparisons = 0;
	for (size_t i = 0; i < spacers.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			string a = spacers[i];
			string b = spacers[j];
			float score = differential_length_string_comparison(a, b);
			score_sum += score;
			comparisons++;
		}
	}


	float mean_score = score_sum / comparisons;

	// divide the score by the number of spacers to punish having more spacers
	return mean_score / spacers.size();

}

void print_spacers(string genome, Crispr crispr, map<string, int> spacer_scores)
{
	
	cout << "\t" << "spacers:" << endl;
	for (string spacer : crispr.spacers)
	{
		printf("\t\t");
		printf("%d/%zd", spacer_scores[spacer], spacer.length());
		// printf(" %d/%zd", spacer_scores[reverse_complement(spacer)], spacer.length()); // unnecessary because BLAST searches both strands.
		printf(" %s\n", spacer.c_str());

	}

	cout << endl;
}


void print_spacers(string genome, Crispr crispr)
{
	map<string, int> spacer_scores;
	for (string spacer : crispr.spacers)
	{
		spacer_scores[spacer] = -1;
	}

	print_spacers(genome, crispr, spacer_scores);

}


void print_repeats(string genome, Crispr crispr, bool reverse_complements)
{
	// string outer = reverse_complements ? "repeats (reverse complements)" : "repeats";
	string outer = "repeats";
	cout << "\t" << outer << endl;

	for (size_t i = 0; i < crispr.repeats.size(); i++)
	{
		// string repeat = reverse_complements ? reverse_complement(repeats[i]) : repeats[i];
		string repeat = crispr.repeats[i];

		int mismatches = mismatch_count(repeat);
		int matches = repeat.length() / 2 - mismatches;
		float score = (float) matches / (float) (repeat.length() / 2);


		int start = crispr.genome_indices[i];
		int end = start + crispr.k - 1;

		if (!POS_STRAND)
		{
			int temp = genome.size() - end;
			end = genome.size() - start;
			start = temp;
		}

		// cout << "\t\t" << matches << "/" << repeat.length() / 2 << " " << score << " " << start << " " << repeat << " " << end << endl;
		printf("\t\t");
		printf("%d/%zd", matches, repeat.length()/2);
		printf(" %d %s %d", start, repeat.c_str(), end);
		printf(" %f\n", score);
	}

	cout << endl;
}

void print_header(string genome, Crispr crispr)
{	
	// vector<string> repeats = get_repeats(crispr, genome);
	// vector<string> spacers = get_spacers(crispr, genome);

	// printf("%d %d %f %f\n", crispr.genome_indices[0], crispr.k, get_conservation_consensus(repeats), get_conservation_spacer(spacers));
	printf("%d %d %f %f\n", crispr.genome_indices[0], crispr.k, crispr.conservation_repeats, crispr.conservation_spacers);
}


void print(string genome, Crispr crispr)
{
	print_header(genome, crispr);
	print_repeats(genome, crispr, false);
	// print_repeats(genome, crispr, true);
	print_spacers(genome, crispr);
}

void print(string genome, Crispr crispr, map<string, int> spacer_scores)
{
	print_header(genome, crispr);
	print_repeats(genome, crispr, false);
	// print_repeats(genome, crispr, true);
	print_spacers(genome, crispr, spacer_scores);
}


void print(string genome, vector<Crispr> crisprs)
{
	printf("printing %zd crisprs\n", crisprs.size());
	for (Crispr crispr : crisprs)
	{
		print(genome, crispr);
	}
}


void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores)
{
	printf("printing %zd crisprs\n", crisprs.size());
	for (Crispr crispr : crisprs)
	{
		print(genome, crispr, spacer_scores);
	}
}




// Is the start and end of the given repeat a subset of any of the repeats of Crispr 'b'? 
bool repeat_substring(Crispr b, int start, int end)
{
	for (size_t i = 0; i < b.genome_indices.size(); i++)
	{
		int repeat_start = b.genome_indices[i];
		int repeat_end = b.genome_indices[i] + b.k - 1;

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

	for (size_t i = 0; i < a.genome_indices.size(); i++)
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

	int a_start = a.genome_indices[0];
	int a_end = a.genome_indices[a.genome_indices.size() - 1] + a.k - 1;

	int b_start = b.genome_indices[0];
	int b_end = b.genome_indices[b.genome_indices.size() - 1] + b.k - 1;



	bool a_before_b = a_start <= b_start;
	bool b_before_a = b_start <= a_start;

	bool a_bleeds_into_b = a_before_b && a_end >= b_start;
	bool b_bleeds_into_a = b_before_a && b_end >= a_start;

	return a_bleeds_into_b || b_bleeds_into_a;
}

bool heuristic_comparison(Crispr a, Crispr b)
{
	// return a.overall_heuristic > b.overall_heuristic;

	if (a.conservation_repeats > b.conservation_repeats)
	{
		return true;
	}
	else if (a.conservation_repeats == b.conservation_repeats)
	{
		return a.conservation_spacers < b.conservation_spacers;
	}
	else
	{
		return false;
	}
	

}