#include "cuda_helpers.h"
#include "cas_helpers.h"
#include "config.h"

#define MAX_MISMATCHES 10

__global__ void
detect_cas_genes(
        ui n_profiles, ui n_profile_kmers, const ui *profile_coords, const ui *profile_kmers,
        ui n_translations, ui n_translation_kmers, const ui *translation_coords, const ui *translation_kmers
) {
    ui profile_coord = blockIdx.x;
    ui translation_coord = threadIdx.x;

    if (profile_coord < n_profiles && translation_coord < n_translations) {
        ui translation_start_index = translation_coords[translation_coord];
        ui translation_end_index =
                translation_coord + 1 >= n_translations ? n_translation_kmers : translation_coords[translation_coord +
                                                                                                   1];
        ui translation_size = translation_end_index - translation_start_index;

        ui profile_start_index = profile_coords[profile_coord];
        ui profile_end_index = profile_coord + 1 >= n_profiles ? n_profile_kmers : profile_coords[profile_coord + 1];
        ui profile_size = profile_end_index - profile_start_index;

        if (profile_size > translation_size) return;

        ui iterations = translation_size - profile_size;

        for (ui iteration = 0; iteration < iterations; iteration++) {
            ui n_mismatch = 0;

            for (ui i = profile_start_index; i < profile_end_index; i++) {
                for (ui j = translation_start_index + iteration;
                     j < translation_start_index + profile_size + iteration; j++) {
                    if (j > n_translation_kmers) break;

                    if (profile_kmers[i] != translation_kmers[j]) {
                        n_mismatch++;
                    }

                    printf("%u : %u\n", translation_kmers[j], profile_kmers[i]);

                    if (n_mismatch > MAX_MISMATCHES) break;
                }

                if (n_mismatch > MAX_MISMATCHES) break;
            }

            if (n_mismatch < MAX_MISMATCHES) {
                printf("Between %u and %u, there were %u mismatches.\n", translation_start_index + iteration,
                       translation_start_index + profile_size + iteration, n_mismatch);
            }
        }
    }
}

vector<Fragment *> Cas::cas_gpu(vector<CasProfile *> &profiles, vector<Translation *> &translations, string &genome) {
    // Steps:
    // 1. get all kmers from each profile (DONE)
    // 2. get all kmers from each translation (DONE)
    // 3. perform comparison on gpu

    ui *h_profile_kmers, *d_profile_kmers;
    ui *h_profile_coords, *d_profile_coords;

    kmer *h_translation_kmers, *d_translation_kmers;
    ui *h_translation_coords, *d_translation_coords;

    AllocationInfo profile_alloc_info = allocate_profiles_on_gpu(
            profiles,
            h_profile_kmers, d_profile_kmers,
            h_profile_coords, d_profile_coords);

    AllocationInfo translation_alloc_info = allocate_translations_on_gpu(
            translations,
            h_translation_kmers, d_translation_kmers,
            h_translation_coords, d_translation_coords);

    ui n_profiles = profile_alloc_info.n_objects;
    ui n_profile_kmers = profile_alloc_info.n_kmers;

    ui n_translations = translation_alloc_info.n_objects;
    ui n_translation_kmers = translation_alloc_info.n_kmers;

    ui numBlocks = n_profiles;
    ui numThreadsPerBlock = n_translations;

    // Perform GPU Kmer comparison operations
    detect_cas_genes<<<numBlocks, numThreadsPerBlock>>>(
            n_profiles, n_profile_kmers, d_profile_coords, d_profile_kmers,
            n_translations, n_translation_kmers, d_translation_coords, d_translation_kmers
    );

    cudaWait();

    // Clean up device memory
    cudaFree(d_profile_kmers);
    cudaFree(d_profile_coords);
    cudaFree(d_translation_kmers);
    cudaFree(d_translation_coords);

    // Clean up host memory
    free(h_profile_kmers);
    free(h_profile_coords);
    free(h_translation_kmers);
    free(h_translation_coords);

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
