#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

#define DEBUG 1

string genome_dir = "T:\\genome";
string domain_table_path = "T:\\truth\\table_domain.tsv";
string serial = "T:\\prospector\\cas\\serial";

std::ofstream out_debug("out_debug.txt");
std::ofstream out_gene("out_gene.txt");
std::ofstream out_domain("out_domain.txt");

void init()
{
    Prospector::device_init();
    CasProfileUtil::load_profiles(serial);
    CasProfileUtil::load_domain_table(domain_table_path);
}

void write_loci(string genome_path, vector<Locus*> loci)
{
    string header = fmt::format("===\t{}\t\t{}\n", genome_path, filesystem::path(genome_path).stem().string());

    out_gene << header;
    out_debug << header;
    out_domain << header;

    for (Locus* l : loci)
    {
        string debug_info = l->to_string_debug();
        string summary_info = l->to_string_summary();

        out_debug << debug_info;

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
}

void prospect_genome(string genome_path)
{
    fmt::print("\n\n\n");

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
    write_loci(genome_path, loci);

    for (Crispr* c : crisprs)
        delete c;

    for (Translation* t : translations)
        delete t;

    for (Fragment* f : fragments)
        delete f;

    start = time(start, "write loci");

    total = time(total, genome_path.c_str());
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
    
    out_debug.close();
    out_gene.close();
    out_domain.close();

    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
