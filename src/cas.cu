#include "cuda_helpers.h"
#include "cas_helpers.h"
#include "config.h"
#include "vector"

__global__ void detect_cas_genes(ui* profile_coords, ui* kmers)
{
    ull i = threadIdx.x + blockDim.x * blockIdx.x;
    printf("This profile starts at: %u\n", profile_coords[i]);
}

::size_t allocate_kmers_on_gpu(vector<CasProfile*>& profiles, ui *&h_kmers, ui *&d_kmers, ui *&h_profile_coords, ui *&d_profile_coords) {
    ::size_t n_profiles = profiles.size();

    ::size_t total_kmers = 0;
    for (const auto &item: profiles)
    {
        total_kmers += item->hash_table.size();
    }

    ull kmers_bytes = total_kmers * sizeof(ui);
    ull profile_locs_bytes = n_profiles * sizeof(ui);

    h_kmers = (ui*) malloc(kmers_bytes);
    h_profile_coords = (ui*) malloc(profile_locs_bytes);

    ui n_kmers = 0;
    for (int i = 0; i < profiles.size(); i++)
    {
        auto hash_table = profiles[i]->hash_table;

        // Store current index of kmer in profile coords array, so we know where this profile starts
        h_profile_coords[i] = n_kmers;

        for (auto const& hash: hash_table)
        {
            h_kmers[n_kmers] = hash;
            n_kmers++;
        }
    }

    // Allocate and copy kmer array
    checkCuda(cudaMalloc(&d_kmers, kmers_bytes));
    checkCuda(cudaMemcpy(d_kmers, h_kmers, kmers_bytes, cudaMemcpyHostToDevice));

    // Allocate and copy profile coords array
    checkCuda(cudaMalloc(&d_profile_coords, profile_locs_bytes));
    checkCuda(cudaMemcpy(d_profile_coords, h_profile_coords, profile_locs_bytes, cudaMemcpyHostToDevice));

    cudaWait();

    return n_profiles;
}

vector<Fragment*> Cas::cas_gpu(vector<CasProfile*>& profiles, vector<Translation*>& translations, string& genome)
{
    // Steps:
    // 1. get all kmers from each profile (DONE)
    // 2. get all kmers from each translation
    // 3. perform comparison on gpu

    ui *h_kmers, *d_kmers;
    ui *h_profile_coords, *d_profile_coords;

    ::size_t n_profiles = allocate_kmers_on_gpu(profiles, h_kmers, d_kmers, h_profile_coords, d_profile_coords);

    detect_cas_genes<<<1, n_profiles>>>(d_profile_coords, d_kmers);

    cudaWait();

    cudaFree(d_kmers);
    cudaFree(d_profile_coords);

    free(h_kmers);
    free(h_profile_coords);

    return {};

//    vector<Fragment*> fragments;
//
//    for (auto profile : profiles)
//    {
//        ui translation_index = -1;
//        for (Translation* t : translations)
//        {
//            translation_index++;
//            vector<ull> index;
//
//            for (ull i = 0; i < t->pure_kmerized_encoded.size(); i++)
//            {
//                bool contains = profile->hash_table.contains(t->pure_kmerized_encoded[i]);
//                if (contains)
//                    index.push_back(i);
//            }
//
//            if (index.size() == 0)
//                continue;
//
//            // dump index to a file for analysis
//            string file_name = Config::path_output / fmt::format("index_dump/{}_{}_{}_{}_{}_{}", CasProfileUtil::domain_table_fetch(profile->identifier), profile->identifier, t->genome_start, t->genome_final, translation_index, index.size());
//
//            if (Config::dump_indices) {
//                std::ofstream out(file_name);
//                for (ull index_location : index)
//                {
//                    out << index_location << "\t" << t->genome_start + index_location << endl;
//                }
//                out.close();
//            }
//
//            vector<vector<ull>> clusters = cluster_index(index);
//            if (!good_clusters(clusters))
//            {
//                continue;
//            }
//            else if (Config::dump_indices)
//            {
//                fmt::print("accepted {}\n", file_name);
//
//                string file_name_clust = fmt::format("{}.clust", file_name);
//                std::ofstream out(file_name_clust);
//                for (vector<ull> clust : clusters)
//                {
//                    for (ull index_location : clust)
//                    {
//                        out << index_location << "\t" << t->genome_start + index_location << endl;
//                    }
//                    out << endl;
//                }
//                out.close();
//            }
//
//
//            Fragment* f = new Fragment;
//            f->reference_genome = &genome;
//            f->reference_crispr = t->reference_crispr;
//            f->reference_translation = t;
//            f->reference_profile = profile;
//            f->clusters = clusters;
//
//            compute_demarc(f);
//
//            if (f->clust_final - f->clust_begin <= 15)
//                continue;
//
//            compute_details(f, genome);
//
//            auto expansion = f->reference_translation->pos ? fragment_expansion_pos : fragment_expansion_neg;
//            expansion(f, genome);
//
//            if (f->genome_begin < 0)
//            {
//                std::exit(1);
//            }
//
//            {
//                fragments.push_back(f);
//            }
//        }
//    }
//
//    vector<Fragment*> filtered_fragments;
//    for (Fragment* f : fragments)
//    {
//        if (CasProfileUtil::domain_table_contains(f->reference_profile->identifier))
//            filtered_fragments.push_back(f);
//        else
//            delete f;
//    }
//
//    std::sort(filtered_fragments.begin(), filtered_fragments.end(), [](Fragment* a, Fragment* b) {return a->expanded_genome_begin < b->expanded_genome_begin; });
//
//    return filtered_fragments;
}
