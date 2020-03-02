#include "util/stdafx.h"
#include "util/util.h"
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




class Profile
{

    public:

        string name;
        string seq;
        vector<string> kmers;

        Profile(string _name, string _path, int _k)
        {
            name = _name;
            seq = parse_fasta(_path).begin()->second;
            kmers = get_kmers(seq, _k); 
        }

        int check_against(vector<string> target_kmers)
        {
            int present = 0;
            for (string query_kmer : kmers)
            {
                for (string target_kmer : target_kmers)
                {
                    int comparison = query_kmer.compare(target_kmer);
                    if (comparison == 0)
                    {
                        present += 1;
                        break;
                    }
                }
            }
            return present;
        }

    private:

};



void cas(string genome, vector<Crispr> crisprs)
{
    clock_t start = clock();

    int k = 5;
    int upstream_size = 10000;    

    vector<Profile> profiles = {
        Profile("thermophilus", "/home/ben/Documents/crispr-data/cas9_amino_thermophilus.fasta", k),
        Profile("pyogenes", "/home/ben/Documents/crispr-data/cas9_amino_pyogenes.fasta", k)
    };

    #pragma omp parallel for      
    for (Crispr crispr : crisprs)
    {
        string upstream = genome.substr(crispr.start - upstream_size, upstream_size);
        vector<string> target_kmers = get_kmers_amino(upstream, k);

        for (Profile profile : profiles)
        {
            int present = profile.check_against(target_kmers);
            printf("profile %s; CRISPR %d %d: %d/%zd\n", profile.name.c_str(), crispr.start, crispr.k, present, profile.kmers.size());

        }
    }



    done(start, "cas detection");

}


vector<Crispr> domain_best(vector<Crispr> crisprs)
{

    printf("filtering %zd crisprs... ", crisprs.size());
    clock_t start = clock();

    
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

map<string, int> get_spacer_scores(vector<Crispr> crisprs)
{
    set<string> all_spacers;
    for (Crispr crispr : crisprs)
    {
        for (string spacer : crispr.spacers)
            all_spacers.insert(spacer);
    }

    map<string, int> spacer_scores = BLAST(all_spacers);
    return spacer_scores;
}



void run()
{        
    map<string, string> genomes = {
        {"thermophilus", parse_fasta("/home/ben/Documents/crispr-data/streptococcus_thermophilus.fasta").begin()->second},
        {"pyogenes", parse_fasta("/home/ben/Documents/crispr-data/pyogenes.fasta").begin()->second}
    };

    string genome = genomes["thermophilus"];

    vector<Crispr> crisprs = prospector_main(genome);
    

    printf("cache crispr information...");
    clock_t cache_start = clock();
    for (size_t i = 0; i < crisprs.size(); i++)
    {
        crisprs[i].update(genome);
    }
    done(cache_start);

    sort(crisprs.begin(), crisprs.end(), heuristic_comparison);

    // debug(genome, crisprs);
    
    vector<Crispr> crisprs_domain_best = domain_best(crisprs);
    
    map<string, int> spacer_scores = get_spacer_scores(crisprs_domain_best);

    vector<Crispr> crisprs_filtered;
    for (Crispr crispr : crisprs_domain_best)
    {
        vector<double> scores;
        for (string spacer : crispr.spacers)
        {
            scores.push_back((double) spacer_scores[spacer] / (double) spacer.size());
        }
        
        if (mean(scores) < 0.5 || crispr.overall_heuristic < 0.5)
        {
            continue;
        }

        crisprs_filtered.push_back(crispr);

    }

    print(genome, crisprs_filtered, spacer_scores);

    cas(genome, crisprs_filtered);
}

int main()
{
    printf("running invoker...\n");
    clock_t start = clock();
    
    run();

    done(start, "invoker");

    return 0;
}

