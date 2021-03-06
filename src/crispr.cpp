
#include "crispr.h"

double get_conservation_consensus(vector<string>& repeats)
{
	string consensus = Util::most_frequent(repeats);

	auto same_length_similarity = [] (string a, string b)
	{
		assert (a.size() == b.size());
		int matches = 0;
		for (ull i = 0; i < a.size(); i++)
			matches += a[i] == b[i] ? 1 : 0;
		return (double) matches / (double) a.length();
	};

	vector<double> similarties;
	similarties.reserve(repeats.size());
	for (ull i = 0; i < repeats.size(); i++)
		similarties.push_back(same_length_similarity(consensus, repeats[i]));
	
	return Util::mean(similarties);
}


double left_aligned_spacer_similarity(string* a, string* b)
{
	ull len = min(a->size(), b->size());
	ull matches = 0;
	for (ull i = 0; i < len; i++)
	{
		matches += a->c_str()[i] == b->c_str()[i] ? 1 : 0;
	}
	return (double) matches / (double) len;
}

double right_aligned_spacer_similarity(string* a, string* b)
{
	string* _big = a;
	string* _small = b;

	if (a->size() < b->size())
	{
		_big = b;
		_small = a;
	}

	ull len = _small->size();
	ull offset = _big->size() - len;
	
	ull matches = 0;
	for (ull i = 0; i < len; i++)
	{
		matches += _big->c_str()[i + offset] == _small->c_str()[i] ? 1 : 0;
	}
	return (double) matches / (double) len;
}

double get_conservation_spacer(vector<string>& spacers)
{
	// compare each spacer against every other space but do not repeat comparisons and do not compare a spacer against itself	
	double score_sum = 0;
	int comparisons = 0;
	for (ull i = 0; i < spacers.size(); i++)
	{
		string* a = &spacers[i];
		for (ull j = 0; j < i; j++)
		{
			string* b = &spacers[j];

			double score = left_aligned_spacer_similarity(a, b);
			if (a->size() != b->size())
			{
				score = (score + right_aligned_spacer_similarity(a, b)) / 2.0;
			}
			score_sum += score;
			comparisons++;
		}
	}
	double mean_score = score_sum / (double) comparisons;
	return mean_score / (double) spacers.size(); // divide the score by the number of spacers to punish having more spacers
}

double get_spacer_variance(vector<string> spacers)
{
	double sum = 0;
	double mean = 0;
	double variance = 0;

	for (ull i = 0; i < spacers.size(); i++)
	{
		sum += (double) spacers[i].length();
	}

	mean = sum / (double) spacers.size();

	for (ull i = 0; i < spacers.size(); i++)
	{
		variance += pow((double) spacers[i].length()    - mean, 2);
	}

	variance = variance / (double) spacers.size();
	return sqrt(variance);
}


// Crispr

Crispr::Crispr(ull k, vector<ull> genome_indices, ull size)
{
	// printf("%d %zd\n", genome_indices[0], size);
	this->k = k;
	this->genome_indices = genome_indices;
	this->size = size;;
}

Crispr::~Crispr()
{
	
}
		
ull Crispr::get_start()
{
	return this->genome_start;
}

ull Crispr::get_final()
{
	return this->genome_final;
}

void Crispr::update(const string& genome)
{
	this->repeats = vector<string>(size);
	this->spacers = vector<string>(size-1);

	for (ull i = 0; i < size; i++)
	{
		repeats[i] = genome.substr(genome_indices[i], k);
	}

	for (ull i = 0; i < size - 1; i++)
	{
		ull current_repeat_end = genome_indices[i] + k;
		ull next_repeat_begin = genome_indices[i+1];
		ull spacer_size = next_repeat_begin - current_repeat_end;
		spacers[i] = genome.substr(current_repeat_end, spacer_size);
	}

	this->genome_start = genome_indices[0];
	this->genome_final = (genome_indices[size-1]) + k - 1;

	this->conservation_repeats = get_conservation_consensus(repeats);
	// this->conservation_spacers = get_conservation_spacer(spacers); 
	this->spacer_variance = get_spacer_variance(spacers) / 100;
	this->conservation_spacers2 = get_conservation_spacer(spacers);

	double subtraction = (50 * conservation_spacers2) + (4*spacer_variance);
	this->overall_heuristic = conservation_repeats - subtraction;
}


//void Crispr::print_generic(const string& genome, function<void(string)>& print_spacer) const
//{
//	// header
//	printf("%d - %d %d\n\n", start, end, k);
//	printf("\t%fh %fr %fs %fv\n\n", overall_heuristic, conservation_repeats, conservation_spacers2, spacer_variance);
//
//	// repeats
//	printf("\trepeats (%zd)\n", repeats.size());
//	printf("\tspacers (%zd)\n", spacers.size());
//
//	for (ull i = 0; i < repeats.size(); i++)
//	{
//		string repeat = repeats[i];
//
//		int mismatches = Util::mismatch_count(repeat);
//		int matches = repeat.length() / 2 - mismatches;
//		double score = (double) matches / (double) (repeat.length() / 2);
//
//		int start = genome_indices[i];
//		int end = start + k - 1;
//
//		int dist = i == 0 ? 0 : genome_indices[i] - (genome_indices[i-1] + k);
//		
//		printf("\t\t");
//		printf("%d/%zd", matches, repeat.length()/2);
//		printf(" %d %s %d", start, repeat.c_str(), end);
//		printf(" %d", dist);
//		printf(" %f", score);
//
//
//
//		if (i < repeats.size() - 1)
//		{
//			string spacer = spacers[i];
//			printf("\t %s\n", spacer.c_str());
//		}
//
//	}
//	cout << endl;
//
//}
//
//void Crispr::print(const string& genome, map<string, int> spacer_scores) const
//{
//	function<void(string)> print_spacer = [&](string spacer) {
//		// printf("%d/%zd", spacer_scores[spacer], spacer.length());
//	};
//
//	print_generic(genome, print_spacer);
//}
//
//void Crispr::print(const string& genome) const
//{
//	function<void(string)> print_spacer = [](string spacer) {
//		// printf("%d/%zd", -1, spacer.length());
//	};
//	print_generic(genome, print_spacer);
//}

string Crispr::identifier_string() const
{
	return fmt::format("{}:{}", this->genome_start, this->k);
}


bool CrisprUtil::heuristic_less(const Crispr* a, const Crispr* b)
{
	return a->overall_heuristic < b->overall_heuristic;
}

bool CrisprUtil::heuristic_greater(const Crispr* a, const Crispr* b)
{
	return a->overall_heuristic > b->overall_heuristic;
}


//void CrisprUtil::print(string genome, vector<Crispr*> crisprs, map<string, int> spacer_scores)
//{
//	printf("printing %zd crisprs\n", crisprs.size());
//	for (Crispr* crispr : crisprs) crispr->print(genome, spacer_scores);
//}
//
//void CrisprUtil::print(string genome, vector<Crispr*> crisprs)
//{
//	printf("printing %zd crisprs\n", crisprs.size());
//	Util::sort(crisprs, CrisprUtil::heuristic_greater);
//	for (Crispr* crispr : crisprs) crispr->print(genome);
//}


// Is the start and end of the given repeat a subset of any of the repeats of Crispr 'b'? 
bool CrisprUtil::repeat_substring(Crispr* b, ull start, ull end)
{
	for (ull i = 0; i < b->size; i++)
	{
		ull repeat_start = b->genome_indices[i];
		ull repeat_end = b->genome_indices[i] + b->k - 1;

		if (start >= repeat_start && end <= repeat_end)
		{
			return true;
		}				
	}
	return false;
}

// Are all the repeats of Crispr 'a' substrings of the repeats in Crispr 'b'?
bool CrisprUtil::repeat_subset(Crispr* a, Crispr* b)
{
	// a cannot be a repeat_subset of b if its k is greater than b
	if (a->k > b->k)
	{
		return false;
	}

	for (ull i = 0; i < a->size; i++)
	{
		if (!repeat_substring(b, a->genome_indices[i], a->genome_indices[i] + a->k - 1))
		{
			return false;
		}
	}

	// all a repeats are substrings of b repeats
	return true;
}

bool CrisprUtil::any_overlap(const Crispr* a, const Crispr* b)
{
	ull a_start = a->genome_indices[0];
	ull a_final = a->genome_indices[a->size - 1] + a->k - 1;
	ull b_start = b->genome_indices[0];
	ull b_final = b->genome_indices[b->size - 1] + b->k - 1;
	return Util::any_overlap(a_start, a_final, b_start, b_final);
}

vector<Crispr*> CrisprUtil::get_domain_best(vector<Crispr*> crisprs)
{
	// this function expects the crisprs to be sorted

    // printf("generating domain best from %zd crisprs... ", crisprs.size());

    // get the best of each domain
    vector<Crispr*> crisprs_domain_best;
    for (ull i = 0; i < crisprs.size(); i++)
    {        
        Crispr* crispr = crisprs[i];

        bool best_already_exists = false;
        for (ull j = 0; j < crisprs_domain_best.size(); j++)
        {
            Crispr* other = crisprs_domain_best[j];
            if (CrisprUtil::any_overlap(crispr, other))
            {
                best_already_exists = true;
                break;
            }
        }

        if (!best_already_exists)
        {
            crisprs_domain_best.push_back(crispr);
        }
    }
    return crisprs_domain_best;
}

vector<Crispr*> CrisprUtil::spacer_score_filtered(vector<Crispr*> crisprs, map<string, int> spacer_scores)
{
    vector<Crispr*> crisprs_filtered;
    for (Crispr* crispr : crisprs)
    {
        vector<double> scores;
        for (string spacer : crispr->spacers)
            scores.push_back((double) spacer_scores[spacer] / (double) spacer.size());

        if (Util::mean(scores) < 0.5)
            continue;

        crisprs_filtered.push_back(crispr);
    }
    return crisprs_filtered;
}

void CrisprUtil::cache_crispr_information(const string& genome, vector<Crispr*>& crisprs)
{
    #pragma omp parallel for
    for (signed long i = 0; i < crisprs.size(); i++)
    {
        crisprs[i]->update(genome);
	}
}



string Crispr::to_string_debug()
{
	std::ostringstream out;

	out << fmt::format("{} - {} {}\n", genome_start, genome_final, k);
	out << fmt::format("{}h {}r {}s {}v\n", overall_heuristic, conservation_repeats, conservation_spacers2, spacer_variance);
	out << fmt::format("\t{} repeats\n", repeats.size());
	out << fmt::format("\t{} spacers\n", spacers.size());

	for (ull i = 0; i < repeats.size(); i++)
	{
		string repeat = repeats[i];
		ull mismatches = Util::mismatch_count(repeat);
		ull matches = repeat.length() / 2 - mismatches;
		double score = (double)matches / (double)(repeat.length() / 2);
		ull start = genome_indices[i];
		ull end = start + k - 1;
		ull dist = i == 0 ? 0 : genome_indices[i] - (genome_indices[i - 1] + k);

		//printf("\t\t");
		//printf("%d/%zd", matches, repeat.length()/2);
		//printf(" %d %s %d", start, repeat.c_str(), end);
		//printf(" %d", dist);
		//printf(" %f", score);

		out << fmt::format("\t\t{} {} {} {}\n", start, repeat, end, i < repeats.size() - 1 ? spacers[i] : "");
	}

	//- 781524	782358	DevR	781329	782373	cas7

	out << "\n";

	return out.str();
}

string Crispr::to_string_summary()
{
	//return fmt::format("{}\t{}\t{}\t{}\t{}h\n", start, end, "?", "CRISPR", overall_heuristic);
	return fmt::format("{}\t{}\t{}\t{}\t{}", genome_start, genome_final, "?", "CRISPR", overall_heuristic);
}