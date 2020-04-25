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

// struct TriFrame
// {
//     ull genome_start;
//     ull genome_end;
//     vector<Translation> translations;
// };

// struct Flanks
// {
//     TriFrame up;
//     TriFrame down;
// };

struct CasProfile
{
    string name;
    string type;
    string raw;
    vector<string> kmers;
    vector<ui> encoded_kmers;
    set<ui> encoded_kmer_set;
    ui* hash_table;
    ui N;
};

struct Fragment
{
    const Crispr* reference_crispr;
    const Translation* reference_translation;
    const CasProfile* reference_profile;
    vector<vector<ull>> clusters;
};

namespace CasUtil
{
    static const string stop = "Z";
    static const char stop_c = 'Z';
    static const ui encoding_size = 5;
    static const ui upstream_size = 10000;
    static const ui k_fragment = 5;

    vector<Translation> get_translations(const string& genome, const vector<Crispr>&);    
    vector<Fragment> cas(const vector<CasProfile>& cas_profiles, const vector<Translation>&);
    void print_fragments(const vector<Crispr>& crisprs, const vector<Fragment>& fragments, const string& genome);
    vector<CasProfile> load(string, ui);
}