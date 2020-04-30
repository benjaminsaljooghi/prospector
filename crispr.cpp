#include "crispr.h"

double get_conservation_consensus(vector<string> repeats)
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
	for (ull i = 0; i < repeats.size(); i++)
		similarties.push_back(same_length_similarity(consensus, repeats[i]));
	
	return Util::mean(similarties);
}


double left_aligned_spacer_similarity(string a, string b)
{
	ull len = min(a.size(), b.size());
	ull matches = 0;
	for (ull i = 0; i < len; i++)
	{
		matches += a[i] == b[i] ? 1 : 0;
	}
	return (double) matches / (double) len;
}

double right_aligned_spacer_similarity(string a, string b)
{
	string _big;
	string _small;

	if (a.size() > b.size())
	{
		_big = a;
		_small = b;
	}
	else
	{
		_big = b;
		_small = a;
	}


	ull len = _small.size();
	ull offset = _big.size() - len;
	
	ull matches = 0;
	for (ull i = 0; i < len; i++)
	{
		matches += _big[i + offset] == _small[i] ? 1 : 0;
	}
	return (double) matches / (double) len;
}

double get_conservation_spacer(vector<string> spacers)
{
	// perform a right aligned and a left aligned spacer mismatch score?
	// if (spacers.size() == 1)
	// {
		// return 0;
	// }

	// compare each spacer against every other space but do not repeat comparisons and do not compare a spacer against itself	
	double score_sum = 0;
	int comparisons = 0;
	for (ull i = 0; i < spacers.size(); i++)
	{
		for (ull j = 0; j < i; j++)
		{
			string a = spacers[i];
			string b = spacers[j];

			double score = left_aligned_spacer_similarity(a, b);
			if (a.size() != b.size())
			{
				score = (score + right_aligned_spacer_similarity(a, b)) / 2.0;
			}
			score_sum += score;
			comparisons++;
		}
	}

	double mean_score = score_sum / (double) comparisons;

	// divide the score by the number of spacers to punish having more spacers
	return mean_score / (double) spacers.size();
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

Crispr::Crispr(unsigned int k, vector<unsigned int> genome_indices, ull size)
{
	// printf("%d %zd\n", genome_indices[0], size);
	this->k = k;
	this->genome_indices = genome_indices;
	this->size = size;;
}
		
void Crispr::update(const string& genome)
{
	this->repeats = vector<string>(size);
	this->spacers = vector<string>(size-1);

	for (ull i = 0; i < size; i++)
	{
		repeats[i] = genome.substr(genome_indices[i], k).c_str();
	}

	for (ull i = 0; i < size - 1; i++)
	{
		unsigned int current_repeat_end = genome_indices[i] + k;
		unsigned int next_repeat_begin = genome_indices[i+1];
		unsigned int spacer_size = next_repeat_begin - current_repeat_end;
		spacers[i] = genome.substr(current_repeat_end, spacer_size);
	}

	this->start = genome_indices[0];
	this->end = (genome_indices[size-1]) + k - 1;

	this->conservation_repeats = get_conservation_consensus(repeats);
	// this->conservation_spacers = get_conservation_spacer(spacers); 
	this->spacer_variance = get_spacer_variance(spacers) / 100;
	this->conservation_spacers2 = get_conservation_spacer(spacers);

	double subtraction = (1.25 * conservation_spacers2) + (4*spacer_variance);
	this->overall_heuristic = conservation_repeats - subtraction;
}


void Crispr::print_generic(const string& genome, function<void(string)>& print_spacer)
{
	// header
	printf("%d - %d %d\n\n", start, end, k);
	printf("\t%fh %fr %fs %fv\n\n", overall_heuristic, conservation_repeats, conservation_spacers2, spacer_variance);

	// repeats
	printf("\trepeats (%zd)\n", repeats.size());

	for (ull i = 0; i < repeats.size(); i++)
	{
		string repeat = repeats[i];

		int mismatches = Util::mismatch_count(repeat);
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
		print_spacer(spacer);
		printf(" %s\n", spacer.c_str());
	}

	cout << endl;

}

void Crispr::print(const string& genome, map<string, int> spacer_scores)
{
	function<void(string)> print_spacer = [&](string spacer) {
		printf("%d/%zd", spacer_scores[spacer], spacer.length());
	};

	print_generic(genome, print_spacer);
}

void Crispr::print(const string& genome)
{
	function<void(string)> print_spacer = [](string spacer) {
		printf("%d/%zd", -1, spacer.length());
	};
	print_generic(genome, print_spacer);
}



// CrisprUtil


bool CrisprUtil::heuristic_less(const Crispr& a, const Crispr& b)
{
	return a.overall_heuristic < b.overall_heuristic;
}

bool CrisprUtil::heuristic_greater(const Crispr& a, const Crispr& b)
{
	return a.overall_heuristic > b.overall_heuristic;
}


void CrisprUtil::print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores)
{
	printf("printing %zd crisprs\n", crisprs.size());
	for (Crispr crispr : crisprs) crispr.print(genome, spacer_scores);
}

void CrisprUtil::print(string genome, vector<Crispr> crisprs)
{
	printf("printing %zd crisprs\n", crisprs.size());
	for (Crispr crispr : crisprs) crispr.print(genome);
}


// Is the start and end of the given repeat a subset of any of the repeats of Crispr 'b'? 
bool CrisprUtil::repeat_substring(Crispr b, unsigned int start, unsigned int end)
{
	for (ull i = 0; i < b.size; i++)
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
bool CrisprUtil::repeat_subset(Crispr a, Crispr b)
{
	// a cannot be a repeat_subset of b if its k is greater than b
	if (a.k > b.k)
	{
		return false;
	}

	for (ull i = 0; i < a.size; i++)
	{
		if (!repeat_substring(b, a.genome_indices[i], a.genome_indices[i] + a.k - 1))
		{
			return false;
		}
	}

	// all a repeats are substrings of b repeats
	return true;
}


bool CrisprUtil::any_overlap(const Crispr& a, const Crispr& b)
{

	ui a_start = a.genome_indices[0];
	ui a_end = a.genome_indices[a.size - 1] + a.k - 1;

	ui b_start = b.genome_indices[0];
	ui b_end = b.genome_indices[b.size - 1] + b.k - 1;


	bool a_before_b = a_start <= b_start;
	bool b_before_a = b_start <= a_start;

	bool a_bleeds_into_b = a_before_b && a_end >= b_start;
	bool b_bleeds_into_a = b_before_a && b_end >= a_start;

	return a_bleeds_into_b || b_bleeds_into_a;
}


vector<Crispr> CrisprUtil::get_domain_best(vector<Crispr> crisprs)
{
	// this function expects the crisprs to be sorted

    // printf("generating domain best from %zd crisprs... ", crisprs.size());
    auto start = time();

    // get the best of each domain
    vector<Crispr> crisprs_domain_best;
    for (ull i = 0; i < crisprs.size(); i++)
    {        
        Crispr crispr = crisprs[i];

        bool best_already_exists = false;
        for (ull j = 0; j < crisprs_domain_best.size(); j++)
        {
            Crispr other = crisprs_domain_best[j];
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
	time(start, "domain best");
    return crisprs_domain_best;
}

vector<Crispr> CrisprUtil::spacer_score_filtered(vector<Crispr> crisprs, map<string, int> spacer_scores)
{
    auto start_time = time();
    vector<Crispr> crisprs_filtered;
    for (Crispr crispr : crisprs)
    {
        vector<double> scores;
        for (string spacer : crispr.spacers)
            scores.push_back((double) spacer_scores[spacer] / (double) spacer.size());

        if (Util::mean(scores) < 0.5)
            continue;

        crisprs_filtered.push_back(crispr);
    }
    time(start_time, "final filtering");
    return crisprs_filtered;
}

void CrisprUtil::cache_crispr_information(const string& genome, vector<Crispr>& crisprs)
{
    auto start = time();
    #pragma omp parallel for
    for (ull i = 0; i < crisprs.size(); i++)
    {
        crisprs[i].update(genome);
	}
    time(start, "cache crisprs");
}


void CrisprUtil::debug(vector<Crispr> crisprs, const string& genome, ui start, ui end)
{

    vector<Crispr> filtered = Util::filter(crisprs, [&](Crispr c) { return c.start > start && c.end < end; } );

    // sort(filtered.begin(), filtered.end(), CrisprUtil::heuristic_less);

    int how_many = filtered.size();
    for (ull i = filtered.size()-how_many; i < filtered.size(); i++)
	{
        filtered[i].print(genome);
	}

	fmt::print("terminating after debug\n");
	exit(0);
}