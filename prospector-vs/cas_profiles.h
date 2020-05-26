

#pragma once

#include <unordered_set>

struct CasProfile
{
    string gn;
    ui N;
    //ui* hash_table;
    unordered_set<ui> hash_table;
    //vector<string> hash_table;

};

namespace CasProfileUtil
{
    static const ui k = 5; // k * encoding size cannot exceed word size.


    //ui get_n(CasProfile& profile);
    //ui gen_n(CasProfile& profile);
    //void load_cache(string);
    //void write_cache(string, vector<CasProfile>);
    //vector<CasProfile> load(string, function<ui(CasProfile&)> get_n);
    vector<CasProfile> cas_profiles_from_uniprot_download(string file_path);
    CasProfile cas_profile_from_tigrfam(string file_path);
    vector<CasProfile> cas_profiles_from_tigrfam(string file_path);
}