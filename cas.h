//std
#include "stdafx.h"
#include <filesystem>
#include "fmt/core.h"
#include "fmt/format.h"
namespace fs = std::filesystem;


//proj
#include "util.h"
#include "crispr.h"



class Translation
{
    public:
        size_t genome_start;
        size_t genome_end;
        map<ui, string> translations_raw;
        map<ui, string> translations_pure;
        map<ui, vector<string>> translations_pure_kmerized;
        map<ui, vector<ui>> translations_pure_kmerized_encoded;
        map<ui, vector<size_t>> pure_mapping;

        Translation(const string& genome, size_t genome_start, size_t genome_end, ui k, bool rc);
        const char* to_string();

        static Translation from_crispr_up(const string& genome, const Crispr& c);
        static Translation from_crispr_down(const string& genonme, const Crispr& c);
};

class CasProfile
{
    public:
        string name;
        string type;
        vector<string> kmers;
        vector<ui> encoded_kmers;
		
		CasProfile(string, ui);

        static vector<CasProfile> load_casprofiles(string, ui);
};

struct Fragment
{
    const Crispr* reference_crispr;
    const Translation* reference_translation;
    const CasProfile* reference_profile;
    vector<vector<size_t>> clusters;
    size_t frame;
};

namespace Cas
{
    vector<Fragment> cas(const string& genome, const vector<Crispr>& crisprs, string cas_dir, const vector<Translation>& downstreams, const vector<Translation>& upstreams);
    void print_fragments(vector<Crispr> crisprs, vector<Fragment> fragments);
}