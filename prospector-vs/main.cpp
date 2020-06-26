#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

#define DEBUG 0

string genome_dir = "T:\\genome";
string domain_table_path = "T:\\prospector-truth\\table_domain.tsv";
string serial = "T:\\prospector\\cas\\serial";
string results_dir = "T:\\results";

void init()
{
    Prospector::device_init();
    CasProfileUtil::load_profiles(serial);
    CasProfileUtil::load_domain_table(domain_table_path);
}

void prospect_genome(string genome_path)
{
    string file_name = filesystem::path(genome_path).stem().string();
    std::filesystem::path results_path = std::filesystem::path(results_dir) / file_name;

    if (std::filesystem::exists(results_path))
    {
        fmt::print("skipping {} because results dir exists\n", genome_path);
        return;
    }

    fmt::print("\n\n");

    auto total = time();
    auto start = time();
    
    string genome = Util::load_genome(genome_path);
    start = time(start, "load genome");

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    fmt::print("returned {} crisprs\n", crisprs.size());
    start = time(start, "get crisprs");

    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    start = time(start, "get translations");

    vector<Fragment*> fragments = Cas::cas(CasProfileUtil::get_profiles(), translations, genome);
    start = time(start, "get fragments");
   
    map<string, bool> crisprs_with_cas;

    for (Crispr* c : crisprs)
        crisprs_with_cas[c->identifier_string()] = false;

    for (Fragment* f : fragments)
    {
        if (f->is_gene())
        {
            crisprs_with_cas[f->reference_crispr->identifier_string()] = true;
        }
    }

    vector<Locus*> loci;
    for (Crispr* c : crisprs)
    {
        if (crisprs_with_cas.at(c->identifier_string()))
        {
            loci.push_back(c);
        }
    }
    for (Fragment* f : fragments) loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) {return a->get_start() < b->get_start(); });

    std::filesystem::create_directory(results_path);

#if DEBUG 1
    std::ofstream out_debug(results_path / "out_debug.txt");
#endif
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_domain(results_path / "out_domain.txt");

    for (Locus* l : loci)
    {


#if DEBUG 1
    string debug_info = l->to_string_debug();
    out_debug << debug_info;
#endif

        string summary_info = l->to_string_summary();


        if (l->is_crispr())
        {
            out_domain << summary_info;
            out_gene << summary_info;
        }

        if (l->is_domain())
        {
            out_domain << summary_info;
        }

        if (l->is_gene())
        {
            out_gene << summary_info;
        }
    }

    start = time(start, "write loci");

    auto total_time = time_diff(total, time());
    total = time(total, "total");

    string finished = fmt::format("// finished in {} ms\n", total_time);

#if DEBUG 1
    out_debug.close();
#endif

    out_gene << finished;
    out_domain << finished;

    out_gene.close();
    out_domain.close();

    for (Crispr* c : crisprs) delete c;
    for (Translation* t : translations) delete t;
    for (Fragment* f : fragments) delete f;
}

void prospect_genome_dir(string genome_dir)
{
    ui i = 0;
    for (const auto& entry : filesystem::directory_iterator(genome_dir))
    {
        if (i++ > 2000) break;
        string genome_path = entry.path().string();
        prospect_genome(genome_path);
    }
}

int main()
{
    auto start_main = time();

    init();

    prospect_genome_dir(genome_dir);
    
    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
