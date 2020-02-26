#include "util/stdafx.h"
#include "util/util.h"

#include "prospector/prospector.h"


map<string, int> BLAST(vector<string> seqs);


int main()
{

    string genome_path("/home/ben/Documents/crispr-data/pyogenes.fasta");
    string genome = parse_fasta(genome_path).begin()->second;

    if (!POS_STRAND)
    {
        genome = reverse_complement(genome);
    }

 
    vector<Crispr> crisprs = prospector_main(genome);



    vector<string> all_spacers;
    for (Crispr crispr : crisprs)
    {
        vector<string> spacers = get_spacers(crispr, genome);
        all_spacers.insert(all_spacers.end(), spacers.begin(), spacers.end());
    }

    map<string, int> spacer_scores = BLAST(all_spacers);


    printf("----------not filtered------------\n");

    print(genome, crisprs, spacer_scores);


    vector<Crispr> crisprs_filtered = vector<Crispr>();
    for (int i = 0; i < crisprs.size(); i++)
    {

        Crispr crispr = crisprs[i];
        vector<string> spacers = get_spacers(crispr, genome);
        vector<string> repeats = get_repeats(crispr, genome);


        // get mean spacer identity

        float spacer_identity_percent_sum = 0;
        for (string spacer : spacers)
        {
            spacer_identity_percent_sum += (float) spacer_scores[spacer] / (float) spacer.length();
        }
        float mean_identity = spacer_identity_percent_sum / spacers.size(); 
    


        // cancel this crispr on the basis of not meeting these requirements

        // insufficient spacer score
        if (mean_identity < 0.5)
        {
            continue;
        }


        // overlaps with another crispr but has a worse conservation than it
        bool overlaps_and_has_worse_conservation = false;
        float this_conservation = repeat_conservation(repeats);
        for (int j = 0; j < crisprs.size(); j++)
        {
            if (i == j)
            {
                continue;
            }

            Crispr other_crispr = crisprs[j];


            if (any_overlap(crispr, other_crispr))
            {
                // these crisprs overlap.

                // cancel the current crispr if it has a worse conservation score than other_crispr.
                // that is, other_crispr is the new candidate bona fide crispr for the current domain.



                float other_conservation = repeat_conservation(get_repeats(other_crispr, genome));

                if (this_conservation < other_conservation)
                {
                    overlaps_and_has_worse_conservation = true;
                    break;
                }

            }
        }
        if (overlaps_and_has_worse_conservation)
        {
            continue;
        }


        // this crispr met all requirements
        crisprs_filtered.push_back(crispr);
    }



    printf("----------filtered----------\n");

    print(genome, crisprs_filtered, spacer_scores);


    return 0;








    // now check upstream of the arrays for cas genes

    int k = 20;
    int upstream_size = 10000;


    // get Cas9
    string cas9_seq = parse_fasta("../crispr-data/addGeneSpCas9.fasta").begin()->second;


    vector<string> cas9_kmers = get_kmers(cas9_seq, k);
    
    string crisrpcasfinder_cas9 = genome.substr(854751, 858857 - 854751);
    vector<string> crisprcasfinder_cas9_kmers = get_kmers(crisrpcasfinder_cas9, k);


    // for (Crispr crispr : crisprs_filtered)
    // {
    //     // transform the upstream 10k into kmers
    //     string upstream = genome.substr(crispr.genome_indices[0] - upstream_size, upstream_size);
    //     vector<string> upstream_kmers = get_kmers(upstream, k);

    //     // where, or do, these kmers overlap? (comparing the set of the cas9_kmers against the set of the upstream kmers)

    //     // what quantity of kmers in the cas9_kmers exists in the upstream_kmers?

    //     int present = 0;

    //     // for (string cas9_kmer : cas9_kmers)
    //     for (string cas9_kmer : crisprcasfinder_cas9_kmers)
    //     {
    //         for (string upstream_kmer : upstream_kmers)
    //         {
    //             present += cas9_kmer.compare(upstream_kmer) == 0 ? 1 : 0;
    //         }
    //     }

    //     printf("for k size %d we have a cas9_kmer count of %zd and a presence of %d\n", crispr.k, crisprcasfinder_cas9_kmers.size(), present);

    // }

}