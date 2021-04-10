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


struct System
{
    vector<Locus*> loci;
    string type;

    ull crispr_count()
    {
        ull count = 0;
        for (Locus* l : loci)
        {
            if (l->is_crispr())
                count++;
        }
        return count;
    }

    ull cas_count()
    {
        ull count = 0;
        for (Locus* l : loci)
        {
            if (!l->is_crispr())
                count++;
        }
        return count;
    }

    void sort_loci()
    {
        std::sort(this->loci.begin(), this->loci.end(), [](Locus* a, Locus* b) { return a->get_start() - b->get_start(); });
    }

    ull get_start()
    {
        Locus* first = this->loci[0];
        return first->get_start();
    }

    ull get_final()
    {
        Locus* last = this->loci[this->loci.size()-1];
        return last->get_final();
    }

    string to_string_summary()
    {
        std::ostringstream out;
        for (Locus* l : this->loci)
            out << l->to_string_summary() << endl;
        return out.str();        
    }

    string to_string_debug()
    {
        std::ostringstream out;
        for (Locus* l : this->loci)
            out << l->to_string_debug() << endl;
        return out.str();        
    }

    bool legitimate_system()
    {
        return this->cas_count() >= 2;
    }



};


vector<MultiFragment*> gen_multifragmetns(vector<Fragment*> fragments)
{
    vector<MultiFragment*> multifragments;
    for (ui i = 0; i < fragments.size(); i++)
    {
        fmt::print("{}: {}\n", i, CasProfileUtil::domain_table_fetch(fragments[i]->reference_profile->identifier));

        fmt::print("multifragment {}\n", i);
        MultiFragment* multifragment = new MultiFragment;
        multifragment->fragments.push_back(fragments[i]);

        for (ui j = i + 1; j < fragments.size(); j++)
        {
            //fmt::print("multifragment comparison {}\n", j);
            //if (fragments[i]->expanded_genome_begin == fragments[j]->expanded_genome_begin &&
                //fragments[i]->expanded_genome_final == fragments[j]->expanded_genome_final)

            bool any_overlap = Util::any_overlap(fragments[i]->genome_begin, fragments[i]->genome_final,
                                                fragments[j]->genome_begin, fragments[j]->genome_final);

            string first = CasProfileUtil::domain_table_fetch(fragments[i]->reference_profile->identifier);
            string second = CasProfileUtil::domain_table_fetch(fragments[j]->reference_profile->identifier);

            bool domain_overlap = (first.find(second) != string::npos) || (second.find(first) != string::npos);

            if (any_overlap && domain_overlap)
            {
                multifragment->fragments.push_back(fragments[j]);
                i = j;
            }
        }

        multifragments.push_back(multifragment);
    }
    return multifragments;
}

vector<System*> gen_systems(vector<Locus*> loci)
{
    vector<System*> systems;

    if (loci.size() == 0)
        return systems;

    System* current = new System;
    current->loci.push_back(loci[0]);
    for (size_t i = 1; i < loci.size(); i++)
    {
        if (loci[i]->get_start() < current->get_final() + 30000)
        {
            current->loci.push_back(loci[i]);
        }
        else
        {
            systems.push_back(current);
            current = new System;
            current->loci.push_back(loci[i]);
        }
    }

    systems.push_back(current);

    vector<System*> filtered_systems;

    for (System* s : systems)
    {
        if (s->legitimate_system())
        {
            filtered_systems.push_back(s);
        }
        else
        {
            delete s;
        }
    }

    return filtered_systems;
}


void prospect_genome(vector<CasProfile*>& profiles, std::filesystem::path genome_path)
{
    fmt::print("\n\n");

    std::filesystem::path results_path = Config::results_dir / genome_path.stem();    
    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_gene_debug(results_path / "out_gene_debug.txt");
    // if (std::filesystem::exists(results_path))
    // {
        // fmt::print("skipping {} because results dir exists\n", genome_path.string());
        // return;
    // }

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Config::crispr_proximal_search ? Cas::crispr_proximal_translations(genome, crisprs) : Cas::get_sixframe(genome, 0, genome.length()-1);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);

    out_gene_debug << fmt::format("BEGIN raw fragments\n");
    for (Fragment* f : fragments)
    {
        out_gene_debug << f->to_string_debug() << endl;   
    }
    out_gene_debug << fmt::format("END raw fragments\n");

    vector<MultiFragment*> multifragments = gen_multifragmetns(fragments);
  
    std::vector<Locus*> loci;

    for (Crispr* c : crisprs)
        loci.push_back(c);

    for (MultiFragment* f : multifragments)
        loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) { return a->get_start() < b->get_start(); });

    vector<System*> systems = gen_systems(loci);

    for (System* system : systems)
    {
        out_gene << system->to_string_summary();
        out_gene_debug << system->to_string_debug() << endl;
    }    

    out_gene.close();
    out_gene_debug.close();

    for (Crispr* c : crisprs) delete c;
    for (Translation* t : translations) delete t;
    for (MultiFragment* m : multifragments) delete m; // this should also delete all fragments
    for (System* s : systems) delete s;
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
    // vector<CasProfile*> filtered = Debug::cas_filter(profiles, "cas9");
   


    // $	542882	545600	-	cas9	cd09643	544647	545600	-1	Cas9_0_II	Cas9_0_II	


    Prospector::device_init();

    // int count = 10;
    // int i = 0;

    // 3057327	3057892	-	cas2,cas2	COG1343,PF09707	

    // profiles = Debug::cas_filter(profiles, "cas2");

    // unordered_set<string> interest{ "GCA_000186245.1_ASM18624v1_genomic.fna" };
    for (const auto& entry : std::filesystem::directory_iterator(Config::genome_dir))
    {
        string filename = entry.path().filename().string();
        // if (interest.contains(filename))
        {
            // Debug::cas_detect(entry, 1578909, 1579230, true, profiles);
            prospect_genome(profiles, entry);
        }
            // Debug::visualize_map(entry.path().string());
            // Debug::cas_detect(entry.path().string(), 3057327, 3057892, false, profiles[0]);        
        
        
        // if (i++ > count)
            // break;
    }

    for (CasProfile* p : profiles)
    {
        delete p;
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

    CasProfileUtil::serialize();
    run();
    // Util::load_genome("/home/ben/crispr/data/genome/assembly/GCF_002863885.1_ASM286388v1_genomic.fna");
    // Debug::cartograph_interpreter(Config::cartograph_prosp, Config::genome_dir);

    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
