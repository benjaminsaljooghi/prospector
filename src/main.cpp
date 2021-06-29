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


vector<MultiFragment*> gen_multifragments(vector<Fragment*> fragments)
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
        // if (true)
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


bool hammer()


bool is_this_seq_what_it_says_it_is(Fragment* fragment) {

    // perform a hmm lookup against the reference profile

    // get reference profile identifier
    auto reference_profile_identifier = fragment->reference_profile->identifier;

    // perform a hammer against this


}


// hardcoded against cas1 for now
vector<ui> hammer(vector<Fragment*>& fragments) {
    
    // runtime:
    // - write candidate sequences to file
    // - system call of hmmscan
    // - parse results
    // - make decisions based on results

    vector<string> seqs;
    for (Fragment* fragment : fragments) {
        seqs.push_back(fragment->protein_sequence);
    } 

    string fasta = Util::seqs_to_fasta(seqs);

    // write candidate sequence to file
    std::ofstream out("/home/ben/hmm/candidates_cas1.txt");
    out << fasta;
    out.close();

    // system call of hmmscan
    string call = "hmmscan --tblout tbloutty.txt /home/ben/hmm/Cas_Cas1.hmm /home/ben/hmm/candidates_cas1.txt > /dev/null";
    system(call.c_str());

    // parse tblout
    std::ifstream infile("tbloutty.txt");    
    string line;
    vector<ui> acceptable_indices;

    while (std::getline(infile, line))
    {
        // fmt::print("line: {}\n", line);
        if (line.starts_with("#"))
            continue;

        vector<string> parse_result = Util::parse(line, " ");
        string query_name = parse_result[2];
        double score = std::stod(parse_result[5]);
        fmt::print("{} {}\n", query_name, score);
        acceptable_indices.push_back(std::stoi(query_name));
    }

    for (ui ting : acceptable_indices)
    {
        fmt::print("{}\n", ting);
    }

    return acceptable_indices;
}

void prospect_genome(vector<CasProfile*>& profiles, std::filesystem::path genome_path)
{
    auto start_prospect = time();

    fmt::print("\n\n");

    std::filesystem::path results_path = Config::results_dir / genome_path.stem();    
    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_gene_debug(results_path / "out_gene_debug.txt");

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Config::crispr_proximal_search ? Cas::crispr_proximal_translations(genome, crisprs) : Cas::get_sixframe(genome, 0, genome.length()-1);
    vector<Fragment*> raw_fragments = Cas::cas(profiles, translations, genome);


    // perform an analysis here of the fragments (either here or multifragments, not sure yet, but let's focus on proof of concept first)
    // let's focus just on cas1 for now
    vector<ui> acceptable_indices = hammer(raw_fragments);

    fmt::print("fasta enumeration of raw fragments:\n");
    for (ui i = 0; i < raw_fragments.size(); i++) {
        fmt::print("{}\n", i);
    }

    fmt::print("we now have acceptable indices\n");
    for (ui index : acceptable_indices) {
        fmt::print("{}\n", index);
    }

    vector<Fragment*> fragments;
    for (ui i = 0; i < raw_fragments.size(); i++) {
        if (Util::contains(acceptable_indices, i))
            fragments.push_back(raw_fragments[i]);
    }


    out_gene_debug << fmt::format("BEGIN raw fragments\n");
    for (Fragment* f : fragments)
    {
        out_gene_debug << f->to_string_debug() << endl;   
    }
    out_gene_debug << fmt::format("END raw fragments\n");

    vector<MultiFragment*> multifragments = gen_multifragments(fragments);
  
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
        // out_gene_debug << system->to_string_debug() << endl;
    }    

    for (Crispr* c : crisprs) delete c;
    for (Translation* t : translations) delete t;
    for (MultiFragment* m : multifragments) delete m;
    for (System* s : systems) delete s;

    auto timed_prospect = time_diff(start_prospect, time());

    out_gene << fmt::format("// finished in {} ms", timed_prospect);

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

void run(vector<CasProfile*>& profiles)
{
    unordered_set<string> interest{ };
    for (const auto& entry : std::filesystem::directory_iterator(Config::genome_dir))
    {
        string filename = entry.path().filename().string();
        if (interest.empty() || (!interest.empty() && interest.contains(filename)))
        {
            prospect_genome(profiles, entry);
        } 
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

    // Debug::cartograph_interpreter(Config::cartograph_prosp, Config::genome_dir); return 0;

    Prospector::device_init();
    // CasProfileUtil::serialize();
    vector<CasProfile*> profiles = CasProfileUtil::deserialize_profiles(Config::serialization_dir);
    vector<CasProfile*> profiles_filtered = Debug::cas_filter(profiles, "cas1");
    
    run(profiles_filtered);
    for (CasProfile* p : profiles) delete p;



    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
