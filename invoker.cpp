#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"




// ProfileExecution




class ProfileExecution
{
    public:
        CasProfile& cas_profile;
        CrisprProfile& crispr_profile;
        map<string, vector<size_t>> result;
		
		// ProfileExecution(CasProfile&, CrisprProfile&);
		// void print();

        		
        ProfileExecution(CasProfile& _cas_profile, CrisprProfile& _crispr_profile)
        :	cas_profile(_cas_profile),
            crispr_profile(_crispr_profile)
        {
            for (auto const& [label, kmers]: crispr_profile.sixway_kmerized)
            {
                for (size_t i = 0; i < kmers.size(); i++)
                {
                    // is this kmer a cas profile kmer?
                    if (contains(cas_profile.kmers, kmers[i]))
                    {
                        // mark this kmer as being a profile kmer
                        result[label].push_back(i);
                    }
                }	
            }
        }


        void print()
        {

            printf("profile %s; CRISPR %d %d\n", cas_profile.name.c_str(), crispr_profile.crispr.start, crispr_profile.crispr.k);


            // compute a full printout of the result

            for (auto const& [label, containment]: result)
            {
                // have current label
                printf("%s \t %zd / %zd (%zd) \n", label.c_str(), containment.size(), crispr_profile.sixway_kmerized[label].size(), cas_profile.kmers.size());


                size_t dist_allowed = 20;
                
                vector<size_t> starts;
                vector<size_t> ends;

                starts.push_back(0);
                for (size_t i = 1; i < containment.size(); i++)
                {
                    size_t prev = containment[i-1];
                    size_t curr = containment[i];

                    size_t dist = curr - prev;
                    if (dist > dist_allowed)
                    {
                        ends.push_back(i-1);
                        starts.push_back(i);
                    }
                }
                ends.push_back(containment.size()-1);

                assert(ends.size() == starts.size());

                for (size_t i = 0; i < starts.size(); i++)
                {
                    size_t start = starts[i];
                    size_t end = ends[i];

                    printf("\t %zd - %zd \t %zd - %zd (%zd)\n", containment[start], containment[end], start, end, end-start+1);

                }



            }
        }

};






void cas(string genome, vector<Crispr> crisprs, const unsigned int k, const size_t upstream_size)
{
    double start = omp_get_wtime();

    vector<CasProfile> cas_profiles = {
        CasProfile("crispr-data/cas9_amino_thermophilus.fasta", k),
        CasProfile("crispr-data/cas9_amino_pyogenes.fasta", k)
    };

    vector<CrisprProfile> crispr_profiles;
    for (Crispr& crispr : crisprs)
    {
        crispr_profiles.push_back(CrisprProfile(genome, crispr, upstream_size, k));
    }

    vector<ProfileExecution> executions;
    #pragma omp parallel for
    for (int i = 0; i < crisprs.size(); i++)
    {
        for (int j = 0; j < cas_profiles.size(); j++)
        {
            ProfileExecution result = ProfileExecution(cas_profiles[j], crispr_profiles[i]);
            #pragma omp critical
            {
                executions.push_back(result);
            }

        }
    }

    for (ProfileExecution& execution : executions)
        execution.print();

    done(start, "cas detection");
}


int main()
{
    printf("running invoker...\n");
    double start = omp_get_wtime();

    map<string, string> genomes =
    {
            {"thermophilus", parse_fasta_single("crispr-data/streptococcus_thermophilus.fasta")},
            {"pyogenes", parse_fasta_single("crispr-data/pyogenes.fasta")}
    };

    string genome = genomes["thermophilus"];

    vector<Crispr> crisprs = Prospector::prospector_main(genome);
    
    CrisprUtil::cache_crispr_information(crisprs, genome);
    
    vector<Crispr> good_heuristic_crisprs;
    for (const Crispr& crispr : crisprs)
    {
        if (crispr.overall_heuristic > 0.5)
            good_heuristic_crisprs.push_back(crispr);
    }

    sort(good_heuristic_crisprs.begin(), good_heuristic_crisprs.end(), CrisprUtil::heuristic_greater);
    
    vector<Crispr> domain_best = CrisprUtil::get_domain_best(good_heuristic_crisprs);

    map<string, int> spacer_scores = CrisprUtil::get_spacer_scores(domain_best);
    
    vector<Crispr> final = CrisprUtil::spacer_score_filtered(domain_best, spacer_scores);

    CrisprUtil::print(genome, final, spacer_scores);

    cas(genome, final, 5, 10000);
    
    done(start, "invoker");

    return 0;
}

