#pragma once

#include "stdafx.h"
#include <filesystem>
#include "fmt/core.h"
#include "fmt/format.h"


#include "util.h"
#include "crispr.h"
#include "cas_profiles.h"


namespace fs = std::filesystem;

struct Translation
{
    const Crispr* reference_crispr;
    ui genome_start;
    ui genome_final;
    bool pos;
    string raw;
    string pure;
    vector<string> pure_kmerized;
    vector<ui> pure_kmerized_encoded;
    vector<ull> pure_mapping;
};

//struct CasSeq
//{
//    string name;
//    string gn;
//    string raw;
//    vector<string> kmers;
//    set<ui> encoded_kmer_set;
//};


struct FragDetails
{
    ull genome_start;
    ull genome_final;
    string genome_translation;  
    size_t quality;
};

struct FragDemarc
{
    ull clust_begin;
    ull clust_final;
};

struct Fragment
{
    const Crispr* reference_crispr;
    const Translation* reference_translation;
    const CasProfile* reference_profile;
    vector<vector<ull>> clusters;

    FragDetails* details;
    FragDemarc* demarc;
};

struct Gene
{
    vector<Fragment> fragments;
    const CasProfile* reference_profile;

    ui size() const
    {
        return fragments[fragments.size()-1].details->genome_final - fragments[0].details->genome_start;
    }
};


namespace Cas
{
    static const ui upstream_size = 10000;
    static const ui cluster_metric_min = 5;
    static const ui max_inter_cluster_dist = 2;

    vector<Translation> get_triframe(const string& genome, ull genome_start, ull genome_final, bool pos);

    vector<Translation> get_sixframe(const string& genome, ull genome_start, ull genome_final);

    vector<Translation> crispr_proximal_translations(const string& genome, const vector<Crispr>&);
    bool* compute_target_map(const vector<const CasProfile*>& cas_profiles, const vector<Translation>& translations);
    vector<Fragment> cas(const vector<const CasProfile*>& cas_profiles, const vector<Translation>&, const string&);
    map<string, vector<Gene>> assemble_genes(const vector<Crispr>& crisprs, const vector<Fragment>& fragments);
    void print_fragment_debug(const Fragment& fragment);
    void print_all(const vector<Crispr>& crisprs, const map<string, vector<Gene>>& crispr_genes, const string& genome);

}

