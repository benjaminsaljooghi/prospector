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
    string raw;
    string pure;
    vector<string> pure_kmerized;
    vector<ui> pure_kmerized_encoded;
    vector<size_t> pure_mapping;
};

struct TriFrame
{
    size_t genome_start;
    size_t genome_end;
    vector<Translation> translations;
};

struct Flanks
{
    TriFrame up;
    TriFrame down;
};

struct CasProfile
{
    string name;
    string type;
    vector<string> kmers;
    vector<ui> encoded_kmers;
};

struct Fragment
{
    const Crispr* reference_crispr;
    const TriFrame* reference_triframe;
    const CasProfile* reference_profile;
    vector<vector<size_t>> clusters;
    size_t frame;
};

namespace CasUtil
{
    vector<Flanks> get_flanks(const string& genome, const vector<Crispr>&);    
    vector<Fragment> cas(const string& genome, const vector<Crispr>& crisprs, const vector<CasProfile>& cas_profiles, const vector<Flanks>&);
    void print_fragments(vector<Crispr> crisprs, vector<Fragment> fragments);
    vector<CasProfile> load(string, ui);
}