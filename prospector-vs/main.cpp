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

//string crispr_type(vector<Gene*>& genes)
//{
//    static const map<string, string> class_lookup
//    {
//        {"I",   "1"},
//        {"III", "1"},
//        {"IV",  "1"},
//        {"II",  "2"},
//        {"V",   "2"},
//        {"VI",  "2"},
//        {"?",   "?"},
//    };
//
//    static const map<string, string> type_lookup
//    {
//        {"cas3",  "I"},
//        {"cas10", "III"},
//        {"cas8",  "IV"},
//        {"cas9",  "II"},
//        {"cas12", "V"},
//        {"cas13", "VI"},
//    };
//
//    for (Gene* gene : genes)
//    {
//        if (type_lookup.contains(gene->gn))
//        {
//            return type_lookup.at(gene->gn);
//        }
//    }
//    return "?";
//}

static vector<CasProfile*> profiles;

void write_all(const string& genome_name, const string& genome, const vector<Crispr*>& crisprs, vector<Fragment*>& fragments, std::ofstream& results)
{
    results << fmt::format("===\t{}\n", genome_name);

    vector<Locus*> loci;
    for (Crispr* c : crisprs) loci.push_back(c);
    for (Fragment* f : fragments) loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) {return a->get_start() < b->get_start(); });


    for (Locus* l : loci)
#if DEBUG == 1
        results << l->to_string_debug();
#else
        results << l->to_string_summary();
#endif

}

void prospect_genome(string genome_path, std::ofstream& results)
{
    fmt::print("\n\n\n");

    auto start = time();

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);

    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    write_all(genome_path, genome, crisprs, fragments, results);

    start = time(start, genome_path.c_str());
}

void prospect_genome_dir(string genome_dir, std::ofstream& results)
{
    ui i = 0;
    for (const auto& entry : filesystem::directory_iterator(genome_dir))
    {
        if (i++ > 3) break;
        string genome_path = entry.path().string();
        prospect_genome(genome_path, results);
    }
}

int main()
{
    // ------ start -----------
    auto start_main = time();


    // ---- serialziation ------
    //CasProfileUtil::pfam_filter("T:\\data\\Pfam-A.seed", "T:\\data\\Pfam-A.filt");
    //vector<const CasProfile*> profiles = CasProfileUtil::pfam("T:\\crispr\\cas\\Pfam-A.full_filt");
    //CasProfileUtil::serialize("T:\\crispr\\cas\\serial_staging", profiles);
    //exit(0);

    // -------- init -----------
    Prospector::device_init();
    profiles = CasProfileUtil::deserialize("T:\\crispr\\cas\\serial");
    std::ofstream results("results.txt");

    // ----------- debug --------
    string debug_path = "T:\\crispr\\supp\\genomes\\GCA_000730285.1_ASM73028v1_genomic.fna";
    
    //Debug::visualize_map(debug_path);
    //prospect_genome(debug_path, results);

    // ----------- prospect --------------
    prospect_genome_dir("T:\\crispr\\supp\\genomes", results);

    // --------- close -----------
    results.close();
    start_main = time(start_main, "main"); return 0;                                                                                                           
}


