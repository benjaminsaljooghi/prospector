//std
#include "stdafx.h"
#include <filesystem>
#include "fmt/core.h"
#include "fmt/format.h"
namespace fs = std::filesystem;


//proj
#include "util.h"
#include "crispr.h"


#define UPSTREAM_SIZE 10000
#define K_FRAGMENT 5

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
    vector<Translation> get_translations(const string& genome, const vector<Crispr>&);    
    vector<Fragment> cas(const vector<CasProfile>& cas_profiles, const vector<Translation>&);
    void print_fragments(const vector<Crispr>& crisprs, const vector<Fragment>& fragments, const string& genome);
    vector<CasProfile> load(string, ui);
}