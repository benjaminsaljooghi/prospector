#include "cuda_helpers.h"
#include "cas_helpers.h"
#include "config.h"
#include "warpcore/warpcore.cuh"
#include "phmap/phmap.h"

#define MAX_MISMATCHES 1000

__global__ void
detect_cas_genes(
        ui n_profiles, ui n_profile_kmers, const ui *profile_coords, const uint64_t *profile_kmers,
        const uint64_t *profile_masks, const ui *profile_masks_coords,
        ui n_translations, ui n_translation_kmers, const ui *translation_coords, const uint64_t *translation_kmers
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

                uint64_t p_k = profile_kmers[profile_start_index + i];
                uint64_t p_mask = profile_masks[profile_start_index + i];
                uint64_t t_k = translation_kmers[translation_start_index + i + iteration];

                uint64_t p_k_masked = p_k & p_mask;
                uint64_t t_k_masked = t_k & p_mask;

                if (p_k_masked != t_k_masked) {
                    n_mismatch++;
                }

                if (n_mismatch > MAX_MISMATCHES) break;
            }

            if (n_mismatch < MAX_MISMATCHES) {
                printf("Between %u and %u, there were %u mismatches.\n", translation_start_index + iteration,
                       translation_start_index + profile_size + iteration, n_mismatch);
                printf("Corresponding Profile Coord - Start %u End %u Size %u\n", profile_start_index,
                       profile_end_index, profile_size);
            }
        }
    }
}


template<typename T>
Cas::AllocationInfo Cas::allocate_kmers_on_gpu(vector<T *> &array,
                                               uint64_t *&h_kmers, uint64_t *&d_kmers,
                                               ui *&h_coords, ui *&d_coords) {
    size_t n = array.size();

    size_t total_kmers = 0;
    for (const auto &item: array) {
        total_kmers += item->binary_kmers.size();
    }

    size_t kmers_bytes = total_kmers * sizeof(uint64_t);
    size_t locs_bytes = n * sizeof(ui);

    h_kmers = (uint64_t *) malloc(kmers_bytes);
    h_coords = (ui *) malloc(locs_bytes);

    printf("Total bytes for array: %lu\n", kmers_bytes + locs_bytes);

    size_t n_kmers = 0;
    for (int i = 0; i < array.size(); i++) {
        // Store current index of kmer in profile coords array, so we know where this profile starts
        h_coords[i] = n_kmers;

        for (auto const &k: array[i]->binary_kmers) {
            h_kmers[n_kmers++] = k;
        }
    }

    // Allocate and copy kmer array
    checkCuda(cudaMalloc(&d_kmers, kmers_bytes));
    checkCuda(cudaMemcpy(d_kmers, h_kmers, kmers_bytes, cudaMemcpyHostToDevice));

    // Allocate and copy profile coords array
    checkCuda(cudaMalloc(&d_coords, locs_bytes));
    checkCuda(cudaMemcpy(d_coords, h_coords, locs_bytes, cudaMemcpyHostToDevice));

    cudaWait();

    AllocationInfo allocationInfo{};
    allocationInfo.n_kmers = n_kmers;
    allocationInfo.n_objects = n;

    return allocationInfo;
}

template<typename T>
Cas::AllocationInfo Cas::allocate_masks_on_gpu(vector<T *> &array,
                                               uint64_t *&h_masks, uint64_t *&d_masks,
                                               ui *&h_mask_coords, ui *&d_mask_coords) {
    size_t n = array.size();

    size_t total_kmers = 0;
    for (const auto &item: array) {
        total_kmers += item->binary_masks.size();
    }

    size_t kmers_bytes = total_kmers * sizeof(uint64_t);
    size_t locs_bytes = n * sizeof(ui);

    h_masks = (uint64_t *) malloc(kmers_bytes);
    h_mask_coords = (ui *) malloc(locs_bytes);

    printf("Total bytes for array: %lu\n", kmers_bytes + locs_bytes);

    size_t n_kmers = 0;
    for (int i = 0; i < array.size(); i++) {
        // Store current index of kmer in profile coords array, so we know where this profile starts
        h_mask_coords[i] = n_kmers;

        for (auto const &k: array[i]->binary_kmers) {
            h_masks[n_kmers++] = k;
        }
    }

    // Allocate and copy kmer array
    checkCuda(cudaMalloc(&d_masks, kmers_bytes));
    checkCuda(cudaMemcpy(d_masks, h_masks, kmers_bytes, cudaMemcpyHostToDevice));

    // Allocate and copy profile coords array
    checkCuda(cudaMalloc(&d_mask_coords, locs_bytes));
    checkCuda(cudaMemcpy(d_mask_coords, h_mask_coords, locs_bytes, cudaMemcpyHostToDevice));

    cudaWait();

    AllocationInfo allocationInfo{};
    allocationInfo.n_kmers = n_kmers;
    allocationInfo.n_objects = n;

    return allocationInfo;
}

vector<Fragment *> Cas::cas_gpu(vector<CasProfile *> &profiles, vector<Translation *> &translations, string &genome) {
    // Steps:
    // 1. get all kmers from each profile (DONE)
    // 2. get all kmers from each translation (DONE)
    // 3. perform comparison on gpu

    uint64_t *h_profile_kmers, *d_profile_kmers;
    ui *h_profile_coords, *d_profile_coords;

    uint64_t *h_profile_masks, *d_profile_masks;
    ui *h_profile_masks_coords, *d_profile_masks_coords;

    uint64_t *h_translation_kmers, *d_translation_kmers;
    ui *h_translation_coords, *d_translation_coords;

    Cas::AllocationInfo profile_alloc_info = Cas::allocate_kmers_on_gpu<CasProfile>(
            profiles,
            h_profile_kmers, d_profile_kmers,
            h_profile_coords, d_profile_coords);

    Cas::allocate_masks_on_gpu<CasProfile>(
            profiles,
            h_profile_masks, d_profile_masks,
            h_profile_masks_coords, d_profile_masks_coords);

    Cas::AllocationInfo translation_alloc_info = Cas::allocate_kmers_on_gpu<Translation>(
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
            n_profiles, n_profile_kmers, d_profile_coords, d_profile_kmers, d_profile_masks, d_profile_masks_coords,
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
