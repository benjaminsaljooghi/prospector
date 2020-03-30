#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"


#define UPSTREAM_SIZE 10000
#define K_FRAGMENT 5


vector<unsigned int> frames{
    0,
    1,
    2
};

// vector<string> pos_labels {
//     "pos_0",
//     "pos_1",
//     "pos_2"
// };

// vector<string> neg_labels {
//     "neg_0",
//     "neg_1",
//     "neg_2"
// };


class ProfileExecution
{
    public:
        CasProfile& cas_profile;
        CrisprProfile& crispr_profile;
        map<unsigned int, vector<size_t>> index;
        map<unsigned int, vector<vector<size_t>>> index_clusters;
        		
        ProfileExecution(CasProfile& _cas_profile, CrisprProfile& _crispr_profile)
        :	cas_profile(_cas_profile),
            crispr_profile(_crispr_profile)
        {
            
        }

        void build_index()
        {
            for (unsigned int frame : frames)
            {
                vector<string> crispr_kmers = this->crispr_profile.translation.translations_pure_kmerized.at(frame);
                vector<string> cas_kmers = this->cas_profile.kmers;

                vector<size_t> indices;
                printf("comparing %zd kmers against %zd kmers\n", crispr_kmers.size(), cas_kmers.size());
                for (size_t i = 0; i < crispr_kmers.size(); i++)
                {
                    if ( contains(cas_kmers, crispr_kmers[i]) )
                    {
                        indices.push_back(i);
                    }
                } 


                if (indices.size() == 0)
                    continue;

                vector<vector<size_t>> clusters; vector<size_t> cluster;

                size_t prev = indices[0];
                for (size_t index : indices)
                {
                    if (index - prev > 5)
                    {
                        vector<size_t> cluster_cp = cluster; clusters.push_back(cluster_cp); cluster.clear();
                    }
                    cluster.push_back(index); prev = index;
                }
                clusters.push_back(cluster);
    


                size_t cluster_requirement = 3;
                bool good_clusters = false;
                for (vector<size_t> cluster : clusters)
                {   
                    if (cluster.size() > cluster_requirement)
                    {
                        good_clusters = true;
                        break;
                    }
                }

                if (good_clusters)
                {
                    index[frame] = indices;
                    index_clusters[frame] = clusters;
                    printf("early termination at frame %d\n", frame);
                    break;
                }

            }

        }

        void single_interpretation(string& genome, unsigned int frame)
        {

            bool contains = index.count(frame) == 1;
            if (!contains)
                return;

            printf("\tframe: %d\n", frame);

            vector<vector<size_t>> clusters = index_clusters[frame];

            // underlying cluster information
            // for (vector<size_t> cluster : clusters)
                // printf("\t \t %zd - %zd (%zd)\n", cluster[0], cluster[cluster.size()-1], cluster.size());

            vector<size_t>  demarc_start;
            vector<size_t>  demarc_end;
            for (vector<size_t> cluster : clusters)
            {
                if (cluster.size() > 1)
                {
                    demarc_start = cluster; break;
                }
            }

            for (size_t i = clusters.size()-1; i >= 0; i--)
            {
                if (clusters[i].size() > 1)
                {
                    demarc_end = clusters[i]; break;
                }
            }

            size_t index_kmer_start = demarc_start[0];
            size_t index_kmer_end = demarc_end[demarc_end.size()-1];

            size_t raw_pos_start = this->crispr_profile.translation.pure_mapping.at(frame)[index_kmer_start];
            size_t raw_pos_end = this->crispr_profile.translation.pure_mapping.at(frame)[index_kmer_end];


            size_t genome_upstream_start = this->crispr_profile.crispr.start - UPSTREAM_SIZE;
            size_t genome_start = genome_upstream_start + (raw_pos_start * 3) + frame;
            size_t genome_end = genome_upstream_start + ((raw_pos_end + K_FRAGMENT) * 3) + frame + 3; // not sure why I need this final +3

            string demarcated_amino = this->crispr_profile.translation.translations_pure.at(frame).substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);
            printf("\t \t %zd -  %zd (%zd - %zd) \n", index_kmer_start, index_kmer_end, genome_start, genome_end );
            printf("%s\n", demarcated_amino.c_str());
            
        }

        void interpret(string genome)
        {

            printf("profile %s; CRISPR %d %d\n", cas_profile.name.c_str(), crispr_profile.crispr.start, crispr_profile.crispr.k);

            for (unsigned int frame : frames) // only iterating labels here
            {
                single_interpretation(genome, frame);
            }

        }

};






void cas(vector<CrisprProfile> crispr_profiles, vector<CasProfile> cas_profiles, string genome)
{
    double start = omp_get_wtime();

    vector<ProfileExecution> executions;
    for (int j = 0; j < cas_profiles.size(); j++)
    {
        for (int i = 0; i < crispr_profiles.size(); i++)
        {
            ProfileExecution result = ProfileExecution(cas_profiles[j], crispr_profiles[i]);
            executions.push_back(result);
        }
    }

    double start_index = omp_get_wtime();
    #pragma omp parallel for
    for (size_t i = 0; i < executions.size(); i++)
    {
        executions[i].build_index();
    }

    
    done(start_index, "build index");

    for (size_t i = 0; i < executions.size(); i++)
    {
        executions[i].interpret(genome);
    }


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
    vector<Crispr> final = CrisprUtil::get_domain_best(good_heuristic_crisprs);
    // map<string, int> spacer_scores = CrisprUtil::get_spacer_scores(domain_best, target_db_path);
    // vector<Crispr> final = CrisprUtil::spacer_score_filtered(domain_best, spacer_scores);
    // CrisprUtil::print(genome, final, spacer_scores);

    CrisprUtil::print(genome, final);

    // cas
    vector<CasProfile> cas_profiles = load_casprofiles(cas_dir, K_FRAGMENT);
    vector<CrisprProfile> crispr_profiles;
    for (Crispr& crispr : final)
    {
        string region = "";
        region += genome.substr(crispr.start - UPSTREAM_SIZE, UPSTREAM_SIZE); 

    // string rc = reverse_complement(nucleotide_sequence);
        // region += genome.substr(crispr.end, UPSTREAM_SIZE); // will need to think about how I'm going to do this. Everything will need to be reversed (all calculations in the interpretation).
        Translation translation(region, K_FRAGMENT);
        CrisprProfile crispr_profile(crispr, translation);
        crispr_profiles.push_back(crispr_profile);
    }
    cas(crispr_profiles, cas_profiles, genome);
    
    done(start, "invoker");

    finish();
    return 0;                                                                                                           
}

