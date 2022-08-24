#include "cuda_include.h"
#include "cas_helpers.h"
#include "config.h"

vector<Fragment*> Cas::cas_gpu(vector<CasProfile*>& profiles, vector<Translation*>& translations, string& genome)
{
    vector<Fragment*> fragments;

    for (auto profile : profiles)
    {
        ui translation_index = -1;
        for (Translation* t : translations)
        {
            translation_index++;
            vector<ull> index;

            for (ull i = 0; i < t->pure_kmerized_encoded.size(); i++)
            {
                bool contains = profile->hash_table.contains(t->pure_kmerized_encoded[i]);
                if (contains)
                    index.push_back(i);
            }

            if (index.size() == 0)
                continue;

            // dump index to a file for analysis
            string file_name = Config::path_output / fmt::format("index_dump/{}_{}_{}_{}_{}_{}", CasProfileUtil::domain_table_fetch(profile->identifier), profile->identifier, t->genome_start, t->genome_final, translation_index, index.size());

            if (Config::dump_indices) {
                std::ofstream out(file_name);
                for (ull index_location : index)
                {
                    out << index_location << "\t" << t->genome_start + index_location << endl;
                }
                out.close();
            }

            vector<vector<ull>> clusters = cluster_index(index);
            if (!good_clusters(clusters))
            {
                continue;
            }
            else if (Config::dump_indices)
            {
                fmt::print("accepted {}\n", file_name);

                string file_name_clust = fmt::format("{}.clust", file_name);
                std::ofstream out(file_name_clust);
                for (vector<ull> clust : clusters)
                {
                    for (ull index_location : clust)
                    {
                        out << index_location << "\t" << t->genome_start + index_location << endl;
                    }
                    out << endl;
                }
                out.close();
            }


            Fragment* f = new Fragment;
            f->reference_genome = &genome;
            f->reference_crispr = t->reference_crispr;
            f->reference_translation = t;
            f->reference_profile = profile;
            f->clusters = clusters;

            compute_demarc(f);

            if (f->clust_final - f->clust_begin <= 15)
                continue;

            compute_details(f, genome);

            auto expansion = f->reference_translation->pos ? fragment_expansion_pos : fragment_expansion_neg;
            expansion(f, genome);

            if (f->genome_begin < 0)
            {
                std::exit(1);
            }

            {
                fragments.push_back(f);
            }
        }
    }

    vector<Fragment*> filtered_fragments;
    for (Fragment* f : fragments)
    {
        if (CasProfileUtil::domain_table_contains(f->reference_profile->identifier))
            filtered_fragments.push_back(f);
        else
            delete f;
    }

    std::sort(filtered_fragments.begin(), filtered_fragments.end(), [](Fragment* a, Fragment* b) {return a->expanded_genome_begin < b->expanded_genome_begin; });

    return filtered_fragments;
}
