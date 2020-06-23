
#pragma once

#include "stdafx.h"
#include <filesystem>
#include "fmt/core.h"
#include "fmt/format.h"


#include "util.h"
#include "crispr.h"
#include "cas_profiles.h"
#include "debug.h"
#include "locus.h"


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

struct Fragment : public Locus
{
    string* reference_genome;

    const Crispr* reference_crispr;
    const Translation* reference_translation;
    const CasProfile* reference_profile;

    vector<vector<ull>> clusters;
    ull clust_begin;
    ull clust_final;

    ull genome_begin;
    ull genome_final;

    ull expanded_genome_begin;
    ull expanded_genome_final;

    ui get_start() { return this->expanded_genome_begin; }
    string to_string_debug();
    string to_string_summary();

    bool is_crispr() { return false;  }
    bool is_domain();
    bool is_gene();
};



namespace Cas
{
    static const ui upstream_size = 20000;
    static const ui cluster_metric_min = 14;
    static const ui max_inter_cluster_dist = 2;

    vector<Translation*> get_triframe(const string& genome, ull genome_start, ull genome_final, bool pos);

    vector<Translation*> get_sixframe(const string& genome, ull genome_start, ull genome_final);

    vector<Translation*> crispr_proximal_translations(const string& genome, vector<Crispr*>&);
    vector<Fragment*> cas(vector<CasProfile*>& cas_profiles, vector<Translation*>&, string&);

    //map<string, string> domain_to_gn;


}

