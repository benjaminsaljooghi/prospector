

#pragma once

#include <unordered_set>
//#include "robin_set.h"
#include "phmap.h"
#include "phmap_dump.h"
#include "stdafx.h"
#include "util.h"
#include <regex>


struct CasProfile
{
    string gn;
    
    
    phmap::flat_hash_set<ui> hash_table;
    

    //ui N;
    //ui* hash_table;
    //unordered_set<ui> hash_table;
    //tsl::robin_set<ui> hash_table;
    //vector<string> hash_table;
    //unordered_set<string> kmer_set;
};



namespace CasProfileUtil
{
    static const ui k = 6; // k * encoding size cannot exceed word size.



    //vector<const CasProfile*> profiles_from_tigrfam_dir(string dir);

    void pfam_filter(string in, string out);
    vector<const CasProfile*> pfam(string path);

    void serialize(string dir, vector<const CasProfile*> profiles);

    vector<CasProfile*> deserialize(string dir);
    void serialize();

}