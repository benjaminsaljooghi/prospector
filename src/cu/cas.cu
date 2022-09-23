#include "cuda_helpers.h"
#include "cas_helpers.h"
#include "config.h"
#include "warpcore/warpcore.cuh"
#include "phmap/phmap.h"


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

            for (ui i = 0; i < profile_size; i++) {
                if (translation_start_index + i + iteration > translation_end_index) break;

                kmer p_k = profile_kmers[profile_start_index + i];
                kmer t_k = translation_kmers[translation_start_index + i + iteration];

                if (p_k != t_k) {
                    n_mismatch++;
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
}
