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


void Crispr::print_generic(string& genome, function<void(string)>& print_spacer)
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
		print_spacer(spacer);
		printf(" %s\n", spacer.c_str());
	}

	cout << endl;

}

void Crispr::print(string& genome, map<string, int> spacer_scores)
{
	function<void(string)> print_spacer = [&](string spacer) {
		printf("%d/%zd", spacer_scores[spacer], spacer.length());
	};

	print_generic(genome, print_spacer);
}

void Crispr::print(string& genome)
{
	function<void(string)> print_spacer = [](string spacer) {
		printf("%d/%zd", -1, spacer.length());
	};

	print_generic(genome, print_spacer);
}



bool CrisprUtil::heuristic_less(const Crispr& a, const Crispr& b)
{
	return a.overall_heuristic < b.overall_heuristic;
}

bool CrisprUtil::heuristic_greater(const Crispr& a, const Crispr& b)
{
	return a.overall_heuristic > b.overall_heuristic;
}



void Crispr::cache_upstream_kmers(string genome, size_t upstream_size, unsigned int _k)
{
	string upstream = genome.substr(start - upstream_size, upstream_size);
	// this->target_kmers = get_kmers_amino(upstream, _k);
}







Profile::Profile(string _name, string _path, unsigned int _k)
{
	name = _name;
	seq = parse_fasta_single(_path);
	kmers = get_kmers(seq, _k);
}

		
ProfileExecution::ProfileExecution(Profile* _profile, Crispr* _crispr)
{
	// how do we demarcate Cas genes?

	profile = _profile;
	crispr = _crispr;

	for (int query = 0; query < profile->kmers.size(); query++)
	{
		string query_kmer = profile->kmers[query];
		for (int target = 0; target < crispr->target_kmers.size(); target++)
		{
			string target_kmer = crispr->target_kmers[target];
			int comparison = query_kmer.compare(target_kmer);
			if (comparison == 0)
			{
				locations_present[query_kmer].push_back(target);
				ordered_positions.push_back(target);
			}
		}
	}

	sort(ordered_positions.begin(), ordered_positions.end());

	hits = ordered_positions.size();
	hits_possible = (profile->kmers).size();
}


void ProfileExecution::print()
{
	printf("profile %s; CRISPR %d %d", profile->name.c_str(), crispr->start, crispr->k);
	printf(": %zd/%zd\n", hits, hits_possible);
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
bool CrisprUtil::repeat_subset(Crispr a, Crispr b)
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


bool CrisprUtil::any_overlap(Crispr a, Crispr b)
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

void CrisprUtil::cas(string genome, vector<Crispr> crisprs, const unsigned int k, const size_t upstream_size)
{
    double start = omp_get_wtime();

    vector<Profile> profiles = {
        Profile("thermophilus", "crispr-data/cas9_amino_thermophilus.fasta", k),
        Profile("pyogenes", "crispr-data/cas9_amino_pyogenes.fasta", k)
    };

    #pragma omp parallel for default(none) shared(crisprs, genome, upstream_size, k)
    for (size_t i = 0; i < crisprs.size(); i++)
    {
        crisprs[i].cache_upstream_kmers(genome, upstream_size, k);
    }
    vector<ProfileExecution> executions;

    #pragma omp parallel for default(none) shared(crisprs, profiles, executions)
    for (int i = 0; i < crisprs.size(); i++)
    {
        for (int j = 0; j < profiles.size(); j++)
        {
            ProfileExecution result = ProfileExecution(&profiles[j], &crisprs[i]);
            #pragma omp critical
            {
                executions.push_back(result);
            };
        }
    }

    for (ProfileExecution& execution : executions)
        execution.print();

    done(start, "cas detection");
}

vector<Crispr> CrisprUtil::get_domain_best(vector<Crispr> crisprs)
{
	// this function expects the crisprs to be sorted

    printf("filtering %zd crisprs... ", crisprs.size());
    double start = omp_get_wtime();

    // get the best of each domain
    vector<Crispr> crisprs_domain_best;
    for (size_t i = 0; i < crisprs.size(); i++)
    {        
        // if there are no overlaps for this crispr (it is unique), then it is a new crispr_best. We have to iterate through all crisprs to find all unique domains.
        Crispr crispr = crisprs[i];

        // check if the domain exists
        bool best_already_exists = false;
        for (size_t j = 0; j < crisprs_domain_best.size(); j++)
        {
            Crispr other = crisprs_domain_best[j];
            if (CrisprUtil::any_overlap(crispr, other))
            {
                // the best of the domain has already been added
                best_already_exists = true;
                break;
            }
        }

        if (!best_already_exists)
        {
            crisprs_domain_best.push_back(crispr);
        }
    }
    done(start);
    return crisprs_domain_best;
}

vector<Crispr> CrisprUtil::score_filtered(vector<Crispr> crisprs, map<string, int> spacer_scores)
{
    double start_time = omp_get_wtime();
    vector<Crispr> crisprs_filtered;
    for (Crispr crispr : crisprs)
    {
        vector<double> scores;
        for (string spacer : crispr.spacers)
            scores.push_back((double) spacer_scores[spacer] / (double) spacer.size());

        if (mean(scores) < 0.5 || crispr.overall_heuristic < 0.5)
            continue;

        crisprs_filtered.push_back(crispr);
    }
    done(start_time, "final filtering");
    return crisprs_filtered;
}

void CrisprUtil::cache_crispr_information(vector<Crispr>& crisprs, string genome)
{
    double start = omp_get_wtime();
    #pragma omp parallel for
    for (size_t i = 0; i < crisprs.size(); i++)
    {
        crisprs[i].update(genome);
    }
    done(start, "cache crispr information");
}

// void CrisprUtil::debug(string genome, vector<Crispr> crisprs)
// {
//     vector<Crispr> of_interest;
//     for (Crispr crispr : crisprs)
//     {
//         if (crispr.start >= DEBUG_START && crispr.end <= DEBUG_END)
//         {
//             of_interest.push_back(crispr);
//         }

//     }
//     sort(of_interest.begin(), of_interest.end(), greater<Crispr>());
//     exit(0);
// }

map<string, int> CrisprUtil::get_spacer_scores(vector<Crispr>& crisprs)
{
    set<string> all_spacers;
    for (Crispr& crispr : crisprs)
        all_spacers.insert(crispr.spacers.begin(), crispr.spacers.end());
    return BLAST(all_spacers);
}
