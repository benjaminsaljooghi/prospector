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
};

namespace CasProfileUtil
{
    static const ull k = 6; // k * encoding size cannot exceed word size.
    void load_domain_map(std::filesystem::path path);
    bool domain_table_contains(string);
    string domain_table_fetch(string);
    vector<CasProfile*> deserialize_profiles(std::filesystem::path);
    void serialize(std::filesystem::path serialization_dir, std::filesystem::path cog_dir);
}