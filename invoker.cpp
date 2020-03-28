#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"




class ProfileExecution
{
    public:
        CasProfile& cas_profile;
        CrisprProfile& crispr_profile;
        map<string, vector<size_t>> result;
        		
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





        void interpret()
        {

            printf("profile %s; CRISPR %d %d\n", cas_profile.name.c_str(), crispr_profile.crispr.start, crispr_profile.crispr.k);


            // compute a full printout of the result

            for (auto const& [label, containment]: result)
            {
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

                
                size_t cluster_requirement = 3;
                bool small_clusters = true;
                for (size_t i = 0; i < starts.size(); i++)
                {   
                    size_t start = starts[i];
                    size_t end = ends[i];
                    size_t len = end-start+1;

                    if (len > cluster_requirement)
                    {
                        small_clusters = false;
                        break;
                    }
                }

                if (small_clusters)
                {
                    continue;
                }



                printf("\t %s \t %zd / %zd (%zd) \n", label.c_str(), containment.size(), crispr_profile.sixway_kmerized[label].size(), cas_profile.kmers.size());



                // print underlying cluster information

                // for (size_t i = 0; i < starts.size(); i++)
                // {
                //     size_t start = starts[i];
                //     size_t end = ends[i];
                //     size_t len = end-start+1;

                //     printf("\t \t %zd - %zd \t %zd - %zd (%zd)\n", containment[start], containment[end], start, end, len);
                // }



                // generate a continuous stretch for which we interpret the gene to exist.
                // It starts at the first non-singleton cluster and ends at the last non-singleton cluster

                // get first non-singleton cluster

                size_t demarc_start = 0;

                for (size_t i = 0; i < starts.size(); i++)
                {
                    size_t start = starts[i];
                    size_t end = ends[i];
                    size_t len = end-start+1;

                    bool non_singleton = len > 1;

                    if (non_singleton)
                    {
                        demarc_start = i;
                        break;
                    }
                
                }

                size_t demarc_end = 0;

                for (size_t i = starts.size()-1; i >= 0; i--)
                {
                    size_t start = starts[i];
                    size_t end = ends[i];
                    size_t len = end-start+1;

                    bool non_singleton = len > 1;

                    if (non_singleton)
                    {
                        demarc_end = i;
                        break;
                    }
                
                }

                // print demarcation
                
                printf("\t \t %zd -  %zd  \n", containment[starts[demarc_start]], containment[ends[demarc_end]]);
                





            }
        }

};






void cas(vector<CrisprProfile> crispr_profiles, vector<CasProfile> cas_profiles)
{
    double start = omp_get_wtime();

    vector<ProfileExecution> executions;
    #pragma omp parallel for
    for (int j = 0; j < cas_profiles.size(); j++)
    {
        for (int i = 0; i < crispr_profiles.size(); i++)
        {
            ProfileExecution result = ProfileExecution(cas_profiles[j], crispr_profiles[i]);
            #pragma omp critical
            {
                executions.push_back(result);
            }
        }
    }

    for (ProfileExecution& execution : executions)
        execution.interpret();

    done(start, "cas detection");
}




void finish()
{
    exit(0);
}

void debug(vector<Crispr> crisprs, string genome)
{
    int how_many = crisprs.size()/2;

    for (size_t i = 0; i < how_many; i++) crisprs[i].print(genome);

    finish();
}



vector<CasProfile> load_casprofiles(string dir, unsigned int k)
{   
    vector<CasProfile> cas_profiles;
    for (const auto& entry : filesystem::directory_iterator(dir))
        cas_profiles.push_back(CasProfile(entry.path(), k));
    return cas_profiles;
}

vector<string> load_genomes(string dir)
{

    vector<string> genomes;
    for (const auto& entry : filesystem::directory_iterator(dir))
        genomes.push_back(parse_fasta_single(entry.path()));
    return genomes;
}


int main()
{
    printf("running invoker...\n");
    double start = omp_get_wtime();

    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";



    string genome = load_genomes(genome_dir)[0];


    vector<Crispr> crisprs = Prospector::prospector_main(genome);
    
    CrisprUtil::cache_crispr_information(crisprs, genome);
    
    vector<Crispr> good_heuristic_crisprs;
    for (const Crispr& crispr : crisprs)
        if (crispr.overall_heuristic > 0.5)
            good_heuristic_crisprs.push_back(crispr);

    sort(good_heuristic_crisprs.begin(), good_heuristic_crisprs.end(), CrisprUtil::heuristic_greater);
    // debug(good_heuristic_crisprs, genome);
    vector<Crispr> domain_best = CrisprUtil::get_domain_best(good_heuristic_crisprs);
    map<string, int> spacer_scores = CrisprUtil::get_spacer_scores(domain_best, target_db_path);
    vector<Crispr> final = CrisprUtil::spacer_score_filtered(domain_best, spacer_scores);
    CrisprUtil::print(genome, final, spacer_scores);



    // cas
    unsigned int k = 5;
    vector<CasProfile> cas_profiles = load_casprofiles(cas_dir, k);
    vector<CrisprProfile> crispr_profiles;
    for (Crispr& crispr : final)
        crispr_profiles.push_back(CrisprProfile(genome, crispr, 10000, k));
    cas(crispr_profiles, cas_profiles);
    



    done(start, "invoker");

    finish();
    return 0;
}

