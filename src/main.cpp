#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time_util.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

std::filesystem::path domain_map_path = "T:/prospector-util/cas/domain_map.tsv";
std::filesystem::path type_table_path = "T:/prospector-util/cas/typing.tsv";
std::filesystem::path serial = "T:/prospector-util/cas/serial_staging";
std::filesystem::path genome_dir = "T:/prospector-util/genome/assembly/";
std::filesystem::path results_dir = "T:/prospector-util/results/";

std::map<string, std::set<string>> type_table;

void load_type_table(std::filesystem::path path)
{
    //std::ifstream table(path);

    //if (!table)
    //{
    //    throw runtime_error(path.string());
    //}

    //string line;
    //while (getline(table, line))
    //{
    //    auto split = Util::parse(line, "\t");
    //    auto type = split[0];
    //    auto required = split[1];
    //    auto additional_any = split[2];
    //    
    //    std::set<string> genes;

    //    for (string gene : Util::parse(, ","))
    //        genes.insert(gene);

    //    type_table[links] = genes;
    //}
}

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

void assert_file(std::filesystem::path path)
{
    if (!std::filesystem::exists(path))
    {
        fmt::print("path does not exist: {}\n", path.string());
        exit(1);
    }
}

void init()
{
    assert_file(domain_map_path);
    assert_file(type_table_path);
    assert_file(serial);
    assert_file(genome_dir);
    assert_file(results_dir);

    Prospector::device_init();
    CasProfileUtil::load_profiles(serial);
    CasProfileUtil::load_domain_table(domain_map_path);
    load_type_table(type_table_path);
}

void prospect_genome(std::filesystem::path genome_path)
{
    std::filesystem::path file_name = genome_path.stem();
    std::filesystem::path results_path = results_dir / file_name;

    if (std::filesystem::exists(results_path))
    {
        fmt::print("skipping {} because results dir exists\n", genome_path.string());
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
    std::sort(fragments.begin(), fragments.end(), [](Fragment* a, Fragment* b) {return a->get_start() < b->get_start(); });
    start = time(start, "get fragments");

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

    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    //for (System* s : systems)
        //out_gene << s->to_string();

    std::vector<Locus*> loci;

    for (Crispr* c : crisprs)
        loci.push_back(c);

    for (Fragment* f : fragments)
        loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) { return a->get_start() < b->get_start(); });

    for (Locus* l : loci)
    {
        fmt::print("locus: {}\n", l->to_string_summary());
        if (!l->is_domain())
            out_gene << l->to_string_summary() << endl;
    }


    start = time(start, "write loci");

    auto total_time = time_diff(total, time());
    total = time(total, "total");

    string finished = fmt::format("// finished in {} ms\n", total_time);

    out_gene << finished;

    out_gene.close();

    for (Crispr* c : crisprs) delete c;
    for (Translation* t : translations) delete t;
    for (Fragment* f : fragments) delete f;
}

void pfam_filt()
{
    std::filesystem::path known_pfam = "T:/prospector-util/pfam_list.txt";
    std::filesystem::path full = "T:/prospector-util/Pfam-A.full";
    std::filesystem::path filt = "T:/prospector-util/Pfam-A.filt";
    CasProfileUtil::pfam_filter(full, filt, known_pfam);
}


unordered_set<string> interest{ "GCF_002355995.1_ASM235599v1_genomic.fna" };

int main()
{
    auto start_main = time();

    init();

    for (const auto& entry : std::filesystem::directory_iterator(genome_dir))
    {
        string filename = entry.path().filename().string();
        if (interest.contains(filename))
            prospect_genome(entry);
    }

    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
