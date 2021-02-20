#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time_util.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"
#include "path.h"

//map<string, System*> system_map;
//for (Fragment* f : fragments)
//{
//    if (!f->is_gene())
//    {
//        continue;
//    }

//    string id = f->reference_crispr->identifier_string();
//    if (system_map.contains(id))
//    {
//        system_map[id]->fragments.push_back(f);
//        continue;
//    }

//    System* system = new System;
//    system->crispr = f->reference_crispr;
//    system->fragments.push_back(f);
//    system_map[id] = system;
//}

//vector<System*> systems;

//for (map<string, System*>::iterator i = system_map.begin(); i != system_map.end(); i++)
//{
//    systems.push_back(i->second);
//}

//for (System* s : systems)
//{
//    s->type = resolve_type(s->fragments);
//}

//std::sort(systems.begin(), systems.end(), [](System* a, System* b) {return a->get_start() - b->get_start(); });
//for (System* s : systems)
    //out_gene << s->to_string();


//std::map<string, std::set<string>> load_type_table(std::filesystem::path path)
//{
//    std::map<string, std::set<string>> type_table;
//
//    std::ifstream table(path);
//
//    if (!table)
//    {
//        throw runtime_error(path.string());
//    }
//
//    string line;
//    while (getline(table, line))
//    {
//        auto split = Util::parse(line, "\t");
//        auto type = split[0];
//        auto required = split[1];
//        auto additional_any = split[2];
//        
//        std::set<string> genes;
//
//        for (string gene : Util::parse(, ","))
//            genes.insert(gene);
//
//        type_table[links] = genes;
//    }
//}

//string resolve_type(vector<Fragment*> fragments)
//{
//    return "N";
//    for (Fragment* f : fragments)
//    {
//        string r = CasProfileUtil::domain_table_fetch(f->reference_profile->identifier);
//        for (const auto& [a, b] : type_table)
//            if (b.contains(r))
//                return a;
//    }
//    return "N";
//}
//
//struct System
//{
//    Crispr* crispr;
//    vector<Fragment*> fragments;
//    string type;
//
//    ui get_start()
//    {
//        return std::min(crispr->get_start(), fragments[0]->get_start());
//    }
//
//    string to_string()
//    {
//        std::ostringstream out;
//
//        vector<Locus*> loci;
//
//
//        loci.push_back(crispr);
//        for (Fragment* f : fragments)
//            loci.push_back(f);
//
//        std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) { return a->get_start() - b->get_start(); });
//
//        for (Locus* l : loci)
//            out << fmt::format("{}\t{}\n", l->to_string_summary(), this->type);
//
//        return out.str();
//    }
//};


namespace Config
{
    std::filesystem::path data_root = "/home/ben/crispr/data/";
    std::filesystem::path util_root = "/home/ben/crispr/prospector-util/";

    // required big data
    std::filesystem::path serialization_dir = data_root / "profiles/";
    std::filesystem::path genome_dir = data_root / "genome/assembly/";

    // required small data
    std::filesystem::path domain_map_path = util_root / "cas/domain_map.tsv";
    std::filesystem::path type_table_path = util_root / "cas/typing.tsv";
    std::filesystem::path results_dir = util_root / "results_prosp/";

    // debug (dev only)
    std::filesystem::path cartograph_prosp = util_root / "cartograph_prosp.tsv";

    bool crispr_proximal_search = true;
}


vector<Fragment*> collect_fragments_crispr_proximal(string& genome, vector<CasProfile*>& profiles, vector<Crispr*>& crisprs)
{
    fmt::print("collecting fragments crispr proximally...\n");
    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    return fragments;
}

vector<Fragment*> collect_fragments_genome_wide(string& genome, vector<CasProfile*>& profiles)
{
    fmt::print("collecting fragments genome wide...\n");
    vector<Translation*> genome_sixframe = Cas::get_sixframe(genome, 0, genome.length()-1);
    vector<Fragment*> fragments = Cas::cas(profiles, genome_sixframe, genome);
    return fragments;
}


void prospect_genome(vector<CasProfile*>& profiles, std::filesystem::path genome_path)
{
    std::filesystem::path results_path = Config::results_dir / genome_path.stem();
    // if (std::filesystem::exists(results_path))
    // {
        // fmt::print("skipping {} because results dir exists\n", genome_path.string());
        // return;
    // }

    // Debug::sage_interpreter("T:\\prospector-util\\report_align.tsv", Config::genome_dir);
    // exit(0);	


    fmt::print("\n\n");

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);

    vector<Fragment*> fragments = Config::crispr_proximal_search ? collect_fragments_crispr_proximal(genome, profiles, crisprs) : collect_fragments_genome_wide(genome, profiles); 

    std::sort(fragments.begin(), fragments.end(), [](Fragment* a, Fragment* b) {return a->expanded_genome_begin < b->expanded_genome_begin; });

    vector<Fragment*> filtered_fragments;

    for (Fragment* f : fragments)
    {
        if (CasProfileUtil::domain_table_contains(f->reference_profile->identifier))
        {
            filtered_fragments.push_back(f);
        }            
    }

    vector<MultiFragment*> multifragments;
    for (ui i = 0; i < filtered_fragments.size(); i++)
    {
        fmt::print("{}: {}\n", i, CasProfileUtil::domain_table_fetch(filtered_fragments[i]->reference_profile->identifier));

        fmt::print("multifragment {}\n", i);
        MultiFragment* multifragment = new MultiFragment;
        multifragment->fragments.push_back(filtered_fragments[i]);

        for (ui j = i + 1; j < filtered_fragments.size(); j++)
        {
            //fmt::print("multifragment comparison {}\n", j);
            //if (filtered_fragments[i]->expanded_genome_begin == filtered_fragments[j]->expanded_genome_begin &&
                //filtered_fragments[i]->expanded_genome_final == filtered_fragments[j]->expanded_genome_final)

            bool any_overlap = Util::any_overlap(filtered_fragments[i]->genome_begin, filtered_fragments[i]->genome_final,
                                                filtered_fragments[j]->genome_begin, filtered_fragments[j]->genome_final);

            string first = CasProfileUtil::domain_table_fetch(filtered_fragments[i]->reference_profile->identifier);
            string second = CasProfileUtil::domain_table_fetch(filtered_fragments[j]->reference_profile->identifier);

            bool domain_overlap = (first.find(second) != string::npos) || (second.find(first) != string::npos);
            // bool domain_overlap = first == second;

            if (any_overlap && domain_overlap)
            {
                multifragment->fragments.push_back(filtered_fragments[j]);
                i = j;
            }
        }

        multifragments.push_back(multifragment);
    }

    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_gene_debug(results_path / "out_gene_debug.txt");

    std::vector<Locus*> loci;

    for (Crispr* c : crisprs)
        loci.push_back(c);

    for (MultiFragment* f : multifragments)
        loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) { return a->get_start() < b->get_start(); });

    for (Locus* l : loci)
    {
        out_gene << l->to_string_summary() << endl;
    }

    for (Locus* l : loci)
    {
        out_gene_debug << l->to_string_debug() << endl;
    }

    out_gene.close();
    out_gene_debug.close();

}

void assert_file(std::filesystem::path path)
{
    if (!std::filesystem::exists(path))
    {
        fmt::print("path does not exist: {}\n", path.string());
        exit(1);
    }
}


void run()
{
    vector<CasProfile*> profiles = CasProfileUtil::deserialize_profiles(Config::serialization_dir);

    for (size_t i = 0; i < profiles.size(); i++)
        fmt::print("loaded: {}\n", profiles[i]->identifier);

    Prospector::device_init();

    // int count = 200;
    // int i = 0;

    // 3057327	3057892	-	cas2,cas2	COG1343,PF09707	

    // profiles = Debug::cas_filter(profiles, "PF09707");

    unordered_set<string> interest{ "GCF_002863665.1_ASM286366v1_genomic.fna" };
    for (const auto& entry : std::filesystem::directory_iterator(Config::genome_dir))
    {
        string filename = entry.path().filename().string();
        if (interest.contains(filename))
            // Debug::cas_detect(entry.path().string(), 3057327, 3057892, false, profiles[0]);        
            prospect_genome(profiles, entry);
        
        // if (i++ > count)
            // break;
    }
}

int main()
{
    auto start_main = time();
    assert_file(Config::domain_map_path);
    assert_file(Config::type_table_path);
    assert_file(Config::serialization_dir);
    assert_file(Config::genome_dir);
    assert_file(Config::results_dir);
    CasProfileUtil::load_domain_map(Config::domain_map_path);    

    run();
    // CasProfileUtil::serialize();
    // Debug::cartograph_interpreter(Config::cartograph_prosp, Config::genome_dir);

    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
