#pragma once

#include <unordered_set>

#include "stdafx.h"
#include "util.h"
#include <regex>

#include "phmap/phmap.h"
#include "phmap/phmap_dump.h"

struct CasProfile
{
    string identifier;
    phmap::flat_hash_set<kmer> hash_table;
    
    ui length_median;
    double length_mean;
    ui length_min;
    ui length_max;

};



namespace CasProfileUtil
{
    static const ull k = 6; // k * encoding size cannot exceed word size.
    void load_domain_map(std::filesystem::path path);
    bool domain_table_contains(string);
    string domain_table_fetch(string);
    bool domain_contained(string domain);
    std::map<string, string> get_domain_map();
    void print_profiles(vector<CasProfile*> profiles);
    vector<CasProfile*> deserialize_profiles(std::filesystem::path);
    void serialize(std::filesystem::path);
}

