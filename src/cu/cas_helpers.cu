#include "cas_helpers.h"
#include "cuda_helpers.h"

::size_t allocate_translations_on_gpu(
        vector<Translation *> &translations,
        kmer *&h_translation_kmers, kmer *&d_translation_kmers,
        ui *&h_translation_coords, ui *&d_translation_coords) {
    ::size_t n_translations = translations.size();

    ::size_t total_kmers = 0;
    for (const auto &translation: translations) {
        total_kmers += translation->pure_kmerized_encoded.size();
    }

    ull kmers_bytes = total_kmers * sizeof(kmer);
    ull translation_locs_bytes = n_translations * sizeof(ui);

    printf("Total bytes for translations: %llu\n", kmers_bytes + translation_locs_bytes);

    h_translation_kmers = (ui *) malloc(kmers_bytes);
    h_translation_coords = (ui *) malloc(translation_locs_bytes);

    ui n_kmers = 0;
    for (int i = 0; i < translations.size(); i++) {
        h_translation_coords[i] = n_kmers;

        for (kmer &k: translations[i]->pure_kmerized_encoded) {
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

    return n_translations;
}
