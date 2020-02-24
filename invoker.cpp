// #include <stdio.h>


// int main()
// {
//     printf("hi\n");
// }



#include "util/stdafx.h"
#include "util/util.h"

#include "prospector/prospector.h"


map<string, int> BLAST(vector<string> seqs);


int main()
{

    // printf("hi\n");

    Util::Prospection prospection = Prospector::prospector_main();

    vector<Util::Locus> crisprs = prospection.crisprs;
    string genome = prospection.genome;
    
    


    vector<string> all_spacers;
    for (Util::Locus crispr : crisprs)
    {
        vector<string> spacers = Util::spacers(genome, crispr);
        all_spacers.insert(all_spacers.end(), spacers.begin(), spacers.end());
    }


    map<string, int> spacer_scores = BLAST(all_spacers);



    
    vector<Util::Locus> crisprs_filtered = vector<Util::Locus>();

    for (Util::Locus crispr : crisprs)
    {

        vector<string> spacers = Util::spacers(genome, crispr);

        float spacer_identity_percent_sum = 0;
        for (string spacer : spacers)
        {
            spacer_identity_percent_sum += spacer_scores[spacer] / spacer.length();
        }

        float mean_identity = spacer_identity_percent_sum / spacers.size(); 
    
        if (mean_identity >= 0.5)
        {
            crisprs_filtered.push_back(crispr);
        }
    }


    for (Util::Locus crispr : crisprs_filtered)
    {
        vector<string> repeats = Util::repeats(genome, crispr);
        vector<string> spacers = Util::spacers(genome, crispr);

        std::ostringstream string_stream;
    
        string_stream << crispr.genome_indices[0] << " " << crispr.k << endl;

        string_stream << "\t" << "repeats:" << endl;
        for (unsigned int i = 0; i < spacers.size(); i++)
        {
            string_stream << "\t\t" << crispr.genome_indices[i] << " " << repeats[i] << " " << crispr.genome_indices[i] + crispr.k << endl;
        }

        string_stream << endl;
    
        string_stream << "\t" << "spacers:" << endl;

        for (string spacer : spacers)
        {
            string_stream << "\t\t" << spacer_scores[spacer] << "/" << spacer.length() << " " << spacer << endl; 
        }

        printf("%s\n\n", string_stream.str().c_str());

    }


    // now check upstream of the arrays for cas genes

    int k = 20;
    int upstream_size = 10000;


    // get Cas9
    string cas9_seq = Util::parse_fasta("../crispr-data/addGeneSpCas9.fasta").begin()->second;


    vector<string> cas9_kmers = Util::get_kmers(cas9_seq, k);
    


    string crisrpcasfinder_cas9 = genome.substr(854751, 858857 - 854751);
    vector<string> crisprcasfinder_cas9_kmers = Util::get_kmers(crisrpcasfinder_cas9, k);

    for (Util::Locus crispr : crisprs_filtered)
    {
        // transform the upstream 10k into kmers
        string upstream = genome.substr(crispr.genome_indices[0] - upstream_size, upstream_size);
        vector<string> upstream_kmers = Util::get_kmers(upstream, k);

        // where, or do, these kmers overlap? (comparing the set of the cas9_kmers against the set of the upstream kmers)

        // what quantity of kmers in the cas9_kmers exists in the upstream_kmers?

        int present = 0;

        // for (string cas9_kmer : cas9_kmers)
        for (string cas9_kmer : crisprcasfinder_cas9_kmers)
        {
            for (string upstream_kmer : upstream_kmers)
            {
                present += cas9_kmer.compare(upstream_kmer) == 0 ? 1 : 0;
            }
        }

        printf("for k size %d we have a cas9_kmer count of %zd and a presence of %d\n", crispr.k, crisprcasfinder_cas9_kmers.size(), present);

    }

}