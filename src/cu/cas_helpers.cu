#include "cas_helpers.h"
#include "cuda_helpers.h"

AllocationInfo allocate_profiles_on_gpu(
        vector<CasProfile *> &profiles,
        uint64_t *&h_profile_kmers, uint64_t *&d_profile_kmers,
        ui *&h_profile_coords, ui *&d_profile_coords) {
    size_t n_profiles = profiles.size();

    size_t total_kmers = 0;
    for (const auto &profile: profiles) {
        total_kmers += profile->binary_kmers.size();
    }

    ull kmers_bytes = total_kmers * sizeof(ui);
    ull profile_locs_bytes = n_profiles * sizeof(ui);

    h_profile_kmers = (uint64_t *) malloc(kmers_bytes);
    h_profile_coords = (ui *) malloc(profile_locs_bytes);

    printf("Total bytes for profiles: %llu\n", kmers_bytes + profile_locs_bytes);

    ui n_kmers = 0;
    for (int i = 0; i < profiles.size(); i++) {
        // Store current index of kmer in profile coords array, so we know where this profile starts
        h_profile_coords[i] = n_kmers;

        for (auto const &k: profiles[i]->binary_kmers) {
            h_profile_kmers[n_kmers++] = k;
        }
    }

    // Allocate and copy kmer array
    checkCuda(cudaMalloc(&d_profile_kmers, kmers_bytes));
    checkCuda(cudaMemcpy(d_profile_kmers, h_profile_kmers, kmers_bytes, cudaMemcpyHostToDevice));

    // Allocate and copy profile coords array
    checkCuda(cudaMalloc(&d_profile_coords, profile_locs_bytes));
    checkCuda(cudaMemcpy(d_profile_coords, h_profile_coords, profile_locs_bytes, cudaMemcpyHostToDevice));

    cudaWait();

    AllocationInfo allocationInfo{};
    allocationInfo.n_kmers = n_kmers;
    allocationInfo.n_objects = n_profiles;

    return allocationInfo;
}

AllocationInfo allocate_translations_on_gpu(
        vector<Translation *> &translations,
        uint64_t *&h_translation_kmers, uint64_t *&d_translation_kmers,
        ui *&h_translation_coords, ui *&d_translation_coords) {
    ::size_t n_translations = translations.size();

    ::size_t total_kmers = 0;
    for (const auto &translation: translations) {
        total_kmers += translation->binary_kmers.size();
    }

    ull kmers_bytes = total_kmers * sizeof(kmer);
    ull translation_locs_bytes = n_translations * sizeof(ui);

    printf("Total bytes for translations: %llu\n", kmers_bytes + translation_locs_bytes);

    h_translation_kmers = (uint64_t *) malloc(kmers_bytes);
    h_translation_coords = (ui *) malloc(translation_locs_bytes);

    ui n_kmers = 0;
    for (int i = 0; i < translations.size(); i++) {
        h_translation_coords[i] = n_kmers;

        for (auto &k: translations[i]->binary_kmers) {
            h_translation_kmers[n_kmers] = k;
            n_kmers++;
        }
    }

    // Allocate and copy kmer array
    checkCuda(cudaMalloc(&d_translation_kmers, kmers_bytes));
    checkCuda(cudaMemcpy(d_translation_kmers, h_translation_kmers, kmers_bytes, cudaMemcpyHostToDevice));

    // Allocate and copy translation coords array
    checkCuda(cudaMalloc(&d_translation_coords, translation_locs_bytes));
    checkCuda(cudaMemcpy(d_translation_coords, h_translation_coords, translation_locs_bytes, cudaMemcpyHostToDevice));

    AllocationInfo allocationInfo{};
    allocationInfo.n_kmers = n_kmers;
    allocationInfo.n_objects = n_translations;

    return allocationInfo;
}
