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
std::ofstream out_debug("out_debug.txt");
std::ofstream out_results("out_results.txt");

string get_accession(string genome_path)
{
    static const std::regex accession_pattern("GC[AF]_[0-9]+\.[0-9]+", regex_constants::icase);
    std::filesystem::path my_path(genome_path);
    string stem = my_path.stem().string();
    std::cmatch match;
    bool result = std::regex_search(stem.c_str(), match, accession_pattern);
    return match.str();
}

string get_header(string genome_path)
{
    return fmt::format("===\t{}\t\t{}\n", genome_path, get_accession(genome_path));
}

void prospect_genome(string genome_path)
{
    fmt::print("\n\n\n");

    auto start = time();

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);

    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    
    start = time(start, genome_path.c_str());

    vector<Locus*> loci;
    for (Crispr* c : crisprs) loci.push_back(c);
    for (Fragment* f : fragments) loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) {return a->get_start() < b->get_start(); });

    out_results << get_header(genome_path);
    for (Locus* l : loci) out_results << l->to_string_summary();
    
    out_debug << get_header(genome_path);
    for (Locus* l : loci) out_debug << l->to_string_debug();
    
}

void prospect_genome_dir(string genome_dir)
{
    ui i = 0;
    for (const auto& entry : filesystem::directory_iterator(genome_dir))
    {
        if (i++ > 3) break;
        string genome_path = entry.path().string();
        prospect_genome(genome_path);
    }
    //return results;
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
    

    // ----------- debug --------
    string debug_path = "T:\\crispr\\supp\\genomes\\GCA_000730285.1_ASM73028v1_genomic.fna";
    
    //Debug::visualize_map(debug_path);
    //prospect_genome(debug_path, results);

    // ----------- prospect --------------
    prospect_genome_dir("T:\\crispr\\supp\\genomes");
    

    // --------- close -----------
    out_debug.close();
    out_results.close();
    start_main = time(start_main, "main"); return 0;                                                                                                           
}


