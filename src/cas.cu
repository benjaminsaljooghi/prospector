#include "cuda_helpers.h"
#include "cas_helpers.h"
#include "config.h"
#include "vector"

struct CudaCasData {
    ull id;
//    phmap::flat_hash_set<kmer> hash_table;
};

__global__ void detect_cas_genes(CudaCasData *data, ::size_t size)
{
    ull i = threadIdx.x + blockDim.x * blockIdx.x;
    printf("ID: %llu with size %lu\n", i, size);

    printf("Current value of data->id: %llu\n", data[i].id);
    data[i].id = i * i;
    printf("After, value of data->id is: %llu\n", data[i].id);
}


vector<Fragment*> Cas::cas_gpu(vector<CasProfile*>& profiles, vector<Translation*>& translations, string& genome)
{
    ull num_profiles = 5; // profiles.size();
    ull bytes = num_profiles * sizeof(CudaCasData);

//    auto *profile_arr = static_cast<CudaCasData *>(malloc(bytes));
//
//    for (int i = 0; i < num_profiles - 1; i++) {
//        profile_arr[i].id = i;
//    }

    CudaCasData *h_profiles, *d_profiles;
    h_profiles = (CudaCasData*) malloc(bytes);
    h_profiles[0].id = 0;
    h_profiles[1].id = 1;
    h_profiles[2].id = 2;
    h_profiles[3].id = 3;
    h_profiles[4].id = 4;

    checkCuda(cudaMalloc(&d_profiles, bytes));
    checkCuda(cudaMemcpy(d_profiles, h_profiles, bytes, cudaMemcpyHostToDevice));

    detect_cas_genes<<<1, num_profiles>>>(d_profiles, num_profiles);

    checkCuda(cudaMemcpy(h_profiles, d_profiles, bytes, cudaMemcpyDeviceToHost));

    cudaWait();

    cudaFree(d_profiles);

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
