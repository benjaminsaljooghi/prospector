

#pragma once

#include <unordered_set>

struct CasProfile
{
    string gn;
    ui N;
    //ui* hash_table;
    unordered_set<ui> hash_table;
    //vector<string> hash_table;
    unordered_set<string> kmer_set;
};

namespace CasProfileUtil
{
    static const ui k = 6; // k * encoding size cannot exceed word size.



    //vector<const CasProfile*> profiles_from_tigrfam_dir(string dir);

    void pfam_filter(string in, string out);
    vector<const CasProfile*> pfam(string path);

}