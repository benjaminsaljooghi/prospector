//#include "util/stdafx.h"
//#include "util/util.h"
#include "prospector/prospector.h"


map<string, int> BLAST(set<string> seqs);



void debug(string genome, vector<Crispr> crisprs)
{
    vector<Crispr> of_interest;
    for (Crispr crispr : crisprs)
    {
        if (crispr.start >= DEBUG_START && crispr.end <= DEBUG_END)
        {
            of_interest.push_back(crispr);
        }

    }
    sort(of_interest.begin(), of_interest.end(), heuristic_comparison);
    print(genome, vector<Crispr>(of_interest.begin(), of_interest.begin() + 5));
    // print(genome, of_interest);
    exit(0);
}






void cas(string genome, vector<Crispr> crisprs, const unsigned int k, const size_t upstream_size)
{
    double start = omp_get_wtime();


    vector<Profile> profiles = {
        Profile("thermophilus", "/home/ben/Documents/crispr-data/cas9_amino_thermophilus.fasta", k),
        Profile("pyogenes", "/home/ben/Documents/crispr-data/cas9_amino_pyogenes.fasta", k)
    };

    #pragma omp parallel for default(none) shared(crisprs, genome, upstream_size, k)
    for (Crispr& crispr : crisprs)
        crispr.cache_upstream_kmers(genome, upstream_size, k);

    size_t size = crisprs.size() * profiles.size();
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


vector<Crispr> domain_best(vector<Crispr> crisprs)
{
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
            if (any_overlap(crispr, other))
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


vector<Crispr> score_filtered(vector<Crispr> crisprs, map<string, int> spacer_scores)
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


map<string, int> get_spacer_scores(vector<Crispr>& crisprs)
{
    set<string> all_spacers;
    for (Crispr& crispr : crisprs)
        all_spacers.insert(crispr.spacers.begin(), crispr.spacers.end());
    return BLAST(all_spacers);
}


void cache_crispr_information(vector<Crispr>& crisprs, string genome)
{
    double start = omp_get_wtime();
    #pragma omp parallel for
    for (Crispr& crispr : crisprs)
        crispr.update(genome);
    done(start, "cache crispr information");
}


int main()
{
    printf("running invoker...\n");
    double start = omp_get_wtime();

    map<string, string> genomes = {
            {"thermophilus", parse_fasta("/home/ben/Documents/crispr-data/streptococcus_thermophilus.fasta").begin()->second},
            {"pyogenes", parse_fasta("/home/ben/Documents/crispr-data/pyogenes.fasta").begin()->second}
    };
    string genome = genomes["thermophilus"];



    vector<Crispr> crisprs = prospector_main(genome);

    cache_crispr_information(crisprs, genome);

    sort(crisprs.begin(), crisprs.end(), heuristic_comparison);

    vector<Crispr> crisprs_domain_best = domain_best(crisprs);

    map<string, int> spacer_scores = get_spacer_scores(crisprs_domain_best);

    vector<Crispr> crisprs_filtered = score_filtered(crisprs_domain_best, spacer_scores);

    print(genome, crisprs_filtered, spacer_scores);
    cas(genome, crisprs_filtered, 5, 10000);



    done(start, "invoker");
    return 0;
}

