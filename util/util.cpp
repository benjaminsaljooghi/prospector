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
	printf("reading %s... ", file_path.c_str());
	clock_t start = clock();
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

	done(start);
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

vector<string> get_kmers(string seq, unsigned int k)
{
	vector<string> kmers;
	for (size_t i = 0; i < seq.length() - k + 1; i++)
	{
		kmers.push_back(seq.substr(i, k));
	}
	return kmers;
}

vector<string> get_kmers_amino(string dna, unsigned int k)
{       
	vector<string> super_amino_acid_kmers;
	for (string amino_acid_seq : sixwaytranslation(dna))
	{
		for (string __kmer : get_kmers(amino_acid_seq, k))
		{
			super_amino_acid_kmers.push_back(__kmer);
		}
	}
	return super_amino_acid_kmers;
}


vector<string> sixwaytranslation(string seq)
{
	map<string, string> codon_table;

	codon_table["TTT"] = 'F';
	codon_table["TTC"] = 'F';
	codon_table["TTA"] = 'L';
	codon_table["TTG"] = 'L';

	codon_table["CTT"] = 'L';
	codon_table["CTC"] = 'L';
	codon_table["CTA"] = 'L';
	codon_table["CTG"] = 'L';
	
	codon_table["ATT"] = 'I';
	codon_table["ATC"] = 'I';
	codon_table["ATA"] = 'I';
	codon_table["ATG"] = 'M';
	
	codon_table["GTT"] = 'V';
	codon_table["GTC"] = 'V';
	codon_table["GTA"] = 'V';
	codon_table["GTG"] = 'V';



	codon_table["TCT"] = "S";
	codon_table["TCC"] = "S";
	codon_table["TCA"] = "S";
	codon_table["TCG"] = "S";
	
	codon_table["CCT"] = "P";
	codon_table["CCC"] = "P";
	codon_table["CCA"] = "P";
	codon_table["CCG"] = "P";
	
	codon_table["ACT"] = "T";
	codon_table["ACC"] = "T";
	codon_table["ACA"] = "T";
	codon_table["ACG"] = "T";
	
	codon_table["GCT"] = "A";
	codon_table["GCC"] = "A";
	codon_table["GCA"] = "A";
	codon_table["GCG"] = "A";



	codon_table["TAT"] = "Y";
	codon_table["TAC"] = "Y";
	codon_table["TAA"] = "STOP"; //STOP
	codon_table["TAG"] = "STOP"; // STOP

	codon_table["CAT"] = "H";
	codon_table["CAC"] = "H";
	codon_table["CAA"] = "Q";
	codon_table["CAG"] = "Q";
	
	codon_table["AAT"] = "N";
	codon_table["AAC"] = "N";
	codon_table["AAA"] = "K";
	codon_table["AAG"] = "K";
	
	codon_table["GAT"] = "D";
	codon_table["GAC"] = "D";
	codon_table["GAA"] = "E";
	codon_table["GAG"] = "E";



	codon_table["TGT"] = "C";
	codon_table["TGC"] = "C";
	codon_table["TGA"] = "STOP"; //STOP
	codon_table["TGG"] = "W";

	codon_table["CGT"] = "R";
	codon_table["CGC"] = "R";
	codon_table["CGA"] = "R";
	codon_table["CGG"] = "R";

	codon_table["AGT"] = "S";
	codon_table["AGC"] = "S";
	codon_table["AGA"] = "R";
	codon_table["AGG"] = "R";

	codon_table["GGT"] = "G";
	codon_table["GGC"] = "G";
	codon_table["GGA"] = "G";
	codon_table["GGG"] = "G";



	// string rc = reverse_complement(seq);

	vector<string> amino_acid_seqs = vector<string>(3);

	int k = 3;

	for (size_t frame = 0; frame < 3; frame++)
	{
		for (size_t i = frame; i + k < seq.size(); i += k)
		{
			string codon = seq.substr(i, k);
			string amino_acid = codon_table[codon];
			if (amino_acid != "STOP")
			{
				amino_acid_seqs[frame] += amino_acid;
			}

		}

	}



	return amino_acid_seqs;
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


double similarity(string a, string b)
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
	
	
	return (double) matches / (double) a.length();
}


double get_conservation_consensus(vector<string> repeats)
{
	string consensus = most_frequent(repeats);

	double similarity_sum = 0;
	for (size_t i = 0; i < repeats.size(); i++)
	{
		similarity_sum += similarity(consensus, repeats[i]);

	}

	return similarity_sum / (double) repeats.size();
}

double bp_match_score(string a, string b)
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

	int max_possible_matches = big_length;

	return (double) matches / (double) max_possible_matches;

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
			double score = bp_match_score(a, b);
			score_sum += score;
			comparisons++;
		}
	}


	double mean_score = score_sum / (double) comparisons;

	// divide the score by the number of spacers to punish having more spacers
	return mean_score / (double) spacers.size();

}

// double get_conservation_spacer_2(vector<string> spacers)
// {
// 	// we want to avoid the presence of there being kmers common across the spacers.

	

// 	// concatenate all the spacers into one sequence and build a frequency list of kmers
	
// 	string superspacer;

// 	for (string spacer : spacers)
// 	{
// 		superspacer += spacer;
// 	}


// 	map<string, int> frequencies;

// 	for (int k = 3; k < 10; k++)
// 	{
// 		for (size_t i = 0; i < superspacer.length() - k + 1; i++)
// 		{
// 			string kmer = superspacer.substr(i, k);
// 			frequencies[kmer] += 1;
// 		}
// 	}


// 	// what is the highest frequency? does the size of k matter?
// 	int max_frequency = 0;
// 	string max_kmer;
// 	for (auto const& x : frequencies)
// 	{
// 		if (x.second > max_frequency)
// 		{	
// 			max_frequency = x.second;
// 			max_kmer = x.first;
// 		}
// 	}

// 	// what does this max_frequency value mean? how do we determine a score from it?



// 	// first, I am interested in this frequency value as a fraction of the number of spacers.
// 	double fraction_a = ((double)1) /( (double) max_frequency / (double) spacers.size());

// 	// second, I am interested in what the length of k is as a fraction of the mean spacer size. The larger k is relative to the mean spacer size, the WORSE the score is.
// 	// double fraction_b = max_frequency / 

// 	return fraction_a;
// }



void print_spacers(string genome, Crispr crispr, map<string, int> spacer_scores)
{
	printf("\tspacers (%zd)\n", crispr.spacers.size());

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
	
	printf("\trepeats (%zd)\n", crispr.repeats.size());

	for (size_t i = 0; i < crispr.repeats.size(); i++)
	{
		// string repeat = reverse_complements ? reverse_complement(repeats[i]) : repeats[i];
		string repeat = crispr.repeats[i];

		int mismatches = mismatch_count(repeat);
		int matches = repeat.length() / 2 - mismatches;
		double score = (double) matches / (double) (repeat.length() / 2);


		int start = crispr.genome_indices[i];
		int end = start + crispr.k - 1;

		// if (!POS_STRAND)
		// {
		// 	int temp = genome.size() - end;
		// 	end = genome.size() - start;
		// 	start = temp;
		// }

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
	printf("%d - %d %d\n\n", crispr.start, crispr.end, crispr.k);
	printf("\t%fh %fr %fs\n\n", crispr.overall_heuristic, crispr.conservation_repeats, crispr.conservation_spacers);
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
bool repeat_substring(Crispr b, unsigned int start, unsigned int end)
{
	for (size_t i = 0; i < b.genome_indices.size(); i++)
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

	unsigned int a_start = a.genome_indices[0];
	unsigned int a_end = a.genome_indices[a.genome_indices.size() - 1] + a.k - 1;

	unsigned int b_start = b.genome_indices[0];
	unsigned int b_end = b.genome_indices[b.genome_indices.size() - 1] + b.k - 1;



	bool a_before_b = a_start <= b_start;
	bool b_before_a = b_start <= a_start;

	bool a_bleeds_into_b = a_before_b && a_end >= b_start;
	bool b_bleeds_into_a = b_before_a && b_end >= a_start;

	return a_bleeds_into_b || b_bleeds_into_a;
}

bool heuristic_comparison(Crispr a, Crispr b)
{

	// bool better_repeat_conservation = a.conservation_repeats > b.conservation_repeats;
	// bool larger_repeat_size = a.k > b.k;


	// if (larger_repeat_size)
	// {
	// 	return true;
	// }

	// if (better_repeat_conservation)
	// {
	// 	return true;
	// }

	// return false;




	// if (a.conservation_repeats == b.conservation_repeats && a.k != b.k)
	// {
	// 	return a.k > b.k;
	// }

	return a.overall_heuristic > b.overall_heuristic;

	// if (a.conservation_repeats > b.conservation_repeats)
	// {
	// 	return true;
	// }
	// else if (a.conservation_repeats == b.conservation_repeats)
	// {
	// 	return a.conservation_spacers < b.conservation_spacers;
	// }
	// else
	// {
	// 	return false;
	// }
	

}