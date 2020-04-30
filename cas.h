#pragma once


//std
#include "stdafx.h"
#include <filesystem>
#include "fmt/core.h"
#include "fmt/format.h"
namespace fs = std::filesystem;

//proj
#include "util.h"
#include "crispr.h"



struct Translation
{
    const Crispr* reference_crispr;
    ui genome_start;
    ui genome_end;
    bool pos;
    string raw;
    string pure;
    vector<string> pure_kmerized;
    vector<ui> pure_kmerized_encoded;
    vector<ull> pure_mapping;
};

struct CasProfile
{
    string name;
    string gn;
    string raw;
    vector<string> kmers;
    set<ui> encoded_kmer_set;
    ui* hash_table;
    ui N;
};

struct FragDetails
{
    ull genome_start;
    ull genome_final;
    string translation;  
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

    // eventually replace these functions with cached versions?
    ui size() const
    {
        return fragments[fragments.size()-1].details->genome_final - fragments[0].details->genome_start;
    }
};



namespace CasUtil
{

    static const ui upstream_size = 10000;
    static const ui k = 5; // k * encoding size cannot exceed word size.
    static const ui cluster_metric_min = 5;

    vector<Translation> get_translations(const string& genome, const vector<Crispr>&);    
    vector<Fragment> cas(const vector<CasProfile>& cas_profiles, const vector<Translation>&, const string&);
    void print_fragments(const vector<Crispr>& crisprs, const vector<Fragment>& fragments, const string& genome);

    ui get_n(CasProfile& profile);
    ui gen_n(CasProfile& profile);


    void load_cache(string);
    void write_cache(string, vector<CasProfile>);
    vector<CasProfile> load(string, function<ui(CasProfile&)> get_n);

}

