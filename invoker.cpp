#include "util/stdafx.h"
#include "util/util.h"

#include "prospector/prospector.h"


map<string, int> BLAST(set<string> seqs);



void run()
{
    string genome_path("/home/ben/Documents/crispr-data/pyogenes.fasta");
    // string genome_path("/home/ben/Documents/crispr-data/pyogenes.fasta");
    
    string genome = parse_fasta(genome_path).begin()->second;

    if (!POS_STRAND)
    {
        genome = reverse_complement(genome);
    }

 
    vector<Crispr> crisprs = prospector_main(genome);

    // printf("--------all crisprs returned by prosector----------");
    // print(genome, crisprs);
    printf("filtering %zd crisprs...\n", crisprs.size());



    // sort crisprs in descending order of heuristic
    sort(crisprs.begin(), crisprs.end(), heuristic_comparison);


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


    // vector<Crispr> crisprs_filtered = vector<Crispr>();


    // #pragma omp parallel for
    // for (size_t i = 0; i < crisprs.size(); i++)
    // {
    //     if (i % 100 == 0)
    //     {
    //         printf("%zd\n", i);
    //     }

    //     Crispr crispr = crisprs[i];

        // // high spacer conservation?
        // if (crispr.conservation_spacers > 0.3)
        // {
        //     continue;
        // }

        // overlaps with another crispr but has a worse conservation than it?
        // bool overlaps_and_has_worse_conservation = false;
        // float this_conservation = crispr.conservation_repeats;
        // for (size_t j = 0; j < crisprs.size(); j++)
        // {
        //     if (i == j)
        //     {
        //         continue;
        //     }

        //     Crispr other_crispr = crisprs[j];

        //     if (any_overlap(crispr, other_crispr))
        //     {
        //         // these crisprs overlap and compete to be the bona fide crispr of this locus.

        //         float other_conservation = other_crispr.conservation_repeats;
        //         if (this_conservation < other_conservation || crispr.conservation_spacers > other_crispr.conservation_spacers)
        //         {
        //             // this_crispr loses candidacy to be the bona fide crispr for this locus.
        //             overlaps_and_has_worse_conservation = true;
        //             break;
        //         }

        //     }
        // }
        // if (overlaps_and_has_worse_conservation)
        // {
        //     continue;
        // }


        // // this crispr met all requirements
        // crisprs_filtered.push_back(crispr);
    // }



    // // get mean spacer identity

    // float spacer_identity_percent_sum = 0;
    // for (string spacer : spacers)
    // {
    //     spacer_identity_percent_sum += (float) spacer_scores[spacer] / (float) spacer.length();
    // }
    // float mean_identity = spacer_identity_percent_sum / spacers.size(); 


    // insufficient spacer score
    // if (mean_identity < 0.5)
    // {
    //     continue;
    // }



    set<string> all_spacers;
    for (Crispr crispr : crisprs_domain_best)
    {
        // all_spacers.insert(all_spacers.end(), spacers.begin(), spacers.end());
        for (string spacer : crispr.spacers)
        { 
            all_spacers.insert(spacer);
            // all_spacers.push_back(reverse_complement(spacer)); // unnecessary because BLAST searches both strands.
        }
    }


    map<string, int> spacer_scores = BLAST(all_spacers);


    // printf("----------not filtered----------\n");
    // print(genome, crisprs, spacer_scores);


    printf("----------filtered----------\n");
    print(genome, crisprs_domain_best, spacer_scores);


    vector<Crispr> of_interest;
    for (Crispr crispr : crisprs)
    {
        if (crispr.k == 36 && crispr.start > 1824944 && crispr.end < 1827552 && crispr.genome_indices.size() > 35)
        {
            of_interest.push_back(crispr);
        }
    }    




    return;


    string spy_cas9_path("/home/ben/Documents/crispr-data/SpyCas9.fasta");
    
    string cas9_seq = parse_fasta(spy_cas9_path).begin()->second;




    int k = 6;
    int upstream_size = 10000;


    vector<string> cas9_kmers = get_kmers(cas9_seq, k);
    

    for (Crispr crispr : crisprs_domain_best)
    {
        // transform the upstream 10k into kmers
        string upstream = genome.substr(crispr.start - upstream_size, upstream_size); // get upstream AND downstream
        vector<string> upstream_kmers = get_kmers(upstream, k);

        // where, or do, these kmers overlap? (comparing the set of the cas9_kmers against the set of the upstream kmers)

        // what quantity of kmers in the cas9_kmers exists in the upstream_kmers?

        int present = 0;

        for (string cas9_kmer : cas9_kmers)
        {
            for (string upstream_kmer : upstream_kmers)
            {
                present += cas9_kmer.compare(upstream_kmer) == 0 ? 1 : 0;
            }
        }

        printf("for CRISPR with start %d and k %d we have a cas9_kmer count of %zd and a presence of %d\n", crispr.start, crispr.k, cas9_kmers.size(), present);

    }


}

int main()
{
    clock_t start = clock();
    
    run();

    done(start, "invoker");

    return 0;
}


// for CRISPR with start 1825009 and k 36 we have a cas9_kmer count of 4098 and a presence of 149
// for CRISPR with start 1085433 and k 36 we have a cas9_kmer count of 4098 and a presence of 350
// for CRISPR with start 1577785 and k 20 we have a cas9_kmer count of 4098 and a presence of 172
// for CRISPR with start 130731 and k 21 we have a cas9_kmer count of 4098 and a presence of 187
// for CRISPR with start 31489 and k 20 we have a cas9_kmer count of 4098 and a presence of 166
// for CRISPR with start 528759 and k 21 we have a cas9_kmer count of 4098 and a presence of 153