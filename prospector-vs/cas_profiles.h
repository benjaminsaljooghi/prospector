

#pragma once

#include <unordered_set>
#include "phmap.h"
#include "phmap_dump.h"
#include "stdafx.h"
#include "util.h"
#include <regex>


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

    void serialize();
}