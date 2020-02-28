#include "util/stdafx.h"
#include "util/util.h"

#include "prospector/prospector.h"


map<string, int> BLAST(set<string> seqs);

void run()
{
    string genome_path("/home/ben/Documents/crispr-data/streptococcus_thermophilus.fasta");
    // string genome_path("/home/ben/Documents/crispr-data/pyogenes.fasta");
    string genome = parse_fasta(genome_path).begin()->second;

    // int start_a = 1079091;
    // int end_a = 1081322;

    // int start_b = 1081340;
    // int end_b = 1083253;

    // Cas9_0_II	1,079,091	1,081,322	+
    // Cas9_1_II	1,081,340	1,083,253

    // string ThCas9 = genome.substr(start_a, end_b - start_a);

    // printf("ThCas9:\n");
    // printf("%s\n", ThCas9.c_str());

    // vector<string> seqs = sixwaytranslation(ThCas9);

    // for (string __seq : seqs)
    // {
    //     printf("frame:\n%s\n", __seq.c_str());
    // }
    

    
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


    // printf("----------top n----------\n");
    // vector<Crispr> top_n(crisprs.begin(), crisprs.begin() + 20);
    // print(genome, top_n, spacer_scores);



    



    printf("----------domain best----------\n");
    print(genome, crisprs_domain_best, spacer_scores);



    vector<Crispr> of_interest;
    for (Crispr crispr : crisprs)
    {
        if (crispr.start > 860000 && crispr.end < 870000)
        {
            of_interest.push_back(crispr);
        }
    }    


    int k = 5;
    int upstream_size = 10000;


    
    map<string, string> profiles;

    profiles["cas9_amino_thermophilus"] = parse_fasta("/home/ben/Documents/crispr-data/cas9_amino_thermophilus.fasta").begin()->second;
    profiles["cas9_amino_pyogenes"] = parse_fasta("/home/ben/Documents/crispr-data/cas9_amino_pyogenes.fasta").begin()->second;


    
	for (auto const& profile_container : profiles)
	{
        string name = profile_container.first;
        string profile = profile_container.second;

        vector<string> query_kmers = get_kmers(profile, k);
        
        #pragma omp parallel for
        for (Crispr crispr : crisprs_domain_best)
        {
            string upstream = genome.substr(crispr.start - upstream_size, upstream_size);
            vector<string> target_kmers = get_kmers_amino(upstream, k);

            // atomic<int> present(0);
            int present = 0;
            for (string query_kmer : query_kmers)
            {
                for (string target_kmer : target_kmers)
                {
                    int to_add = query_kmer.compare(target_kmer) == 0 ? 1 : 0;
                    present += to_add;
                }
            }

            printf("for profile %s, for CRISPR with start %d and k %d we have a query kmer count of %zd and a presence of %d\n", name.c_str(), crispr.start, crispr.k, query_kmers.size(), present);

        }
    
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