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
    phmap::flat_hash_set<ui> hash_table;
};


namespace CasProfileUtil
{
    static const ui k = 6; // k * encoding size cannot exceed word size.

    void load_profiles(std::filesystem::path);

    vector<CasProfile*>& get_profiles();

    void load_domain_table(std::filesystem::path);

    bool domain_table_contains(string);

    string domain_table_fetch(string);

    void pfam_filter(std::filesystem::path in, std::filesystem::path out, std::filesystem::path known_path);

    vector<const CasProfile*> profiles_from_tigrfam_dir(std::filesystem::path tigr_table, std::filesystem::path tigr_dir);

    void serialize();
}