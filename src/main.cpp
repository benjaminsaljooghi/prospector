#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time_util.h"
#include "cas_profiles.h"
#include "array_discovery.h"
#include "path.h"
#include <boost/algorithm/string.hpp>
#include "config.h"

// // hardcoded against cas1 for now
// vector<ui> hammer(vector<Fragment*> fragments) {
    
//     // runtime:
//     // - write candidate sequences to file
//     // - system call of hmmscan
//     // - parse results
//     // - make decisions based on results

//     vector<string> seqs;
//     for (Fragment* fragment : fragments) {
//         seqs.push_back(fragment->protein_sequence);
//     } 

//     string fasta = Util::seqs_to_fasta(seqs);

//     // write candidate sequence to file
//     std::ofstream out("/home/ben/crispr/hmm/candidates_cas1.txt");
//     out << fasta;
//     out.close();

//     // system call of hmmscan
//     string call = "hmmscan --tblout tbloutty.txt /home/ben/crispr/hmm/Cas_Cas1.hmm /home/ben/crispr/hmm/candidates_cas1.txt > /dev/null";
//     system(call.c_str());

//     // parse tblout
//     std::ifstream infile("tbloutty.txt");    
//     string line;
//     vector<ui> acceptable_indices;

//     while (std::getline(infile, line))
//     {
//         // fmt::print("line: {}\n", line);
//         if (line.starts_with("#"))
//             continue;

//         vector<string> parse_result = Util::parse(line, " ");
//         string query_name = parse_result[2];
//         double score = std::stod(parse_result[5]);
//         fmt::print("{} {}\n", query_name, score);
//         acceptable_indices.push_back(std::stoi(query_name));
//     }

//     for (ui ting : acceptable_indices)
//     {
//         fmt::print("{}\n", ting);
//     }

//     return acceptable_indices;
// }

// vector<Fragment*> acceptable_indices_filter_methodology(vector<Fragment*> raw_fragments)
// {
//     // perform an analysis here of the fragments (either here or multifragments, not sure yet, but let's focus on proof of concept first)
//     // let's focus just on cas1 for now
//     vector<ui> acceptable_indices = hammer(raw_fragments);

//     fmt::print("fasta enumeration of raw fragments:\n");
//     for (ui i = 0; i < raw_fragments.size(); i++) {
//         fmt::print("{}\n", i);
//     }

//     fmt::print("we now have acceptable indices\n");
//     for (ui index : acceptable_indices) {
//         fmt::print("{}\n", index);
//     }

//     vector<Fragment*> fragments;
//     for (ui i = 0; i < raw_fragments.size(); i++) {
//         if (Util::contains(acceptable_indices, i))
//             fragments.push_back(raw_fragments[i]);
//     }

//     return fragments;
// }



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

vector<string> get_all_hmm_files()
{
    vector<string> paths;

    // step 1: iterate over the dir:
    std::filesystem::path directory = "/home/ben/crispr/prospector-data/bin_hmm_makarova";

    // does the split[0] of the stem match the cas_type in a case insensitive way?
    for (const auto& entry : std::filesystem::directory_iterator(directory))
	{
		string filename = entry.path().filename().string();
        vector<string> extensions = Util::parse(filename, ".");
        string last_extension = extensions[extensions.size()-1];
        
        if (last_extension != "hmm")
            continue;

        string cas_type_of_file = Util::parse(filename, "_")[0];
        paths.push_back(entry.path().string());
    }
    return paths;
}

vector<string> get_hmm_files(string cas_type)
{
    vector<string> paths;

    // step 1: iterate over the dir:
    std::filesystem::path directory = "/home/ben/crispr/prospector-data/bin_hmm_makarova";

    // does the split[0] of the stem match the cas_type in a case insensitive way?
    for (const auto& entry : std::filesystem::directory_iterator(directory))
	{
		string filename = entry.path().filename().string();
        vector<string> extensions = Util::parse(filename, ".");
        string last_extension = extensions[extensions.size()-1];
        
        if (last_extension != "hmm")
            continue;

        string cas_type_of_file = Util::parse(filename, "_")[0];
        if (boost::iequals(cas_type_of_file, cas_type))
        {
            paths.push_back(entry.path().string());
        }
    }
    return paths;
}

string get_hmm_file(string hmm_identifier)
{
    fmt::print("seeking HMM file with identifier: {}\n", hmm_identifier);
    return fmt::format("/home/ben/crispr/hmm/hmm_files/{}.HMM", hmm_identifier);
}

bool have_hmm_file(string hmm_file)
{
    bool exists = std::filesystem::exists(hmm_file);
    if (!exists)
    {
        fmt::print("does NOT exist: {}\n", hmm_file);
    }
    return exists;
}

bool singleton_hammer(string& singleton_seq, string hmm_file)
{
    if (!have_hmm_file(hmm_file))
    {
        throw new exception();
    }


    vector<string> singleton_vec;
    singleton_vec.push_back(singleton_seq);
    string fasta = Util::seqs_to_fasta(singleton_vec);

    std::filesystem::path temp_fasta_file = "/media/ben/temp/crispr_lfs/temp_fasta.txt";

    // write candidate sequence to file
    std::ofstream out(temp_fasta_file.string());
    out << fasta;
    out.close();

    string call = fmt::format("hmmscan --cpu 10 --tblout tbloutty.txt {} {} > /dev/null", hmm_file, temp_fasta_file.string());
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

    // return true;
    assert(acceptable_indices.size() <= 1);
    return acceptable_indices.size() > 0;
}

// bool singleton_hammer(Fragment* fragment) {
//     string hmm_file;
    
//     // attempt an exact match
//     hmm_file = get_hmm_file(fragment->reference_profile->identifier);
//     if (!have_hmm_file(hmm_file)) {
//         // exact match could be found, so this is probably like a COG profile. In which case we get the profile's domain and then we find a suitable HMM for that domain.
//     }

//     return singleton_hammer(fragment->protein_sequence, hmm_file);
// }

bool multimatch_hammer(Fragment* fragment) {
    string cas_type_of_fragment = CasProfileUtil::domain_table_fetch(fragment->reference_profile->identifier);
    vector<string> hmm_files = get_hmm_files(cas_type_of_fragment);

    for (string hmm_file : hmm_files) {
        if (singleton_hammer(fragment->protein_sequence, hmm_file)) {
            return true;
        }
    }
    return false;
}

// vector<Fragment*> individuated_singleton_methodology(vector<Fragment*> raw_fragments) {
//     vector<Fragment*> fragments;
//     for (Fragment* fragment : raw_fragments) {
//         if (singleton_hammer(fragment)) {
//             fragments.push_back(fragment);
//         }
//         else {
//             fmt::print("rejecting {} {} {} {} on the basis of hammer\n", fragment->reference_profile->identifier, CasProfileUtil::domain_table_fetch(fragment->reference_profile->identifier), fragment->expanded_genome_begin, fragment->expanded_genome_final);
//         }
//     }
//     return fragments;
// }

vector<Fragment*> multimatch_singleton_methodology(vector<Fragment*> raw_fragments) {
    vector<Fragment*> fragments;
    for (Fragment* fragment : raw_fragments) {
        if (multimatch_hammer(fragment)) {
            fragments.push_back(fragment);
        }
        else {
            fmt::print("rejecting {} {} {} {} on the basis of hammer\n", fragment->reference_profile->identifier, CasProfileUtil::domain_table_fetch(fragment->reference_profile->identifier), fragment->expanded_genome_begin, fragment->expanded_genome_final);
        }
    }
    return fragments;
}

void prospect_genome(vector<CasProfile*>& profiles, std::filesystem::path genome_path)
{
    auto start_prospect = time();

    fmt::print("\n\n");

    std::filesystem::path results_path = Config::path_results / genome_path.stem();    
    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_gene_debug(results_path / "out_gene_debug.txt");

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Config::crispr_proximal_search ? Cas::crispr_proximal_translations(genome, crisprs) : Cas::get_sixframe(genome, 0, genome.length()-1);
    vector<Fragment*> raw_fragments = Cas::cas(profiles, translations, genome);

    out_gene_debug << fmt::format("BEGIN raw fragments\n");
    for (Fragment* f : raw_fragments)
    {
        out_gene_debug << f->to_string_debug() << endl;   
    }
    out_gene_debug << fmt::format("END raw fragments\n");

//    vector<Fragment*> fragments = multimatch_singleton_methodology(raw_fragments);
//    vector<MultiFragment*> multifragments = gen_multifragments(fragments);
  
    std::vector<Locus*> loci;

    for (Crispr* c : crisprs)
        loci.push_back(c);

//    for (MultiFragment* f : multifragments)
//        loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) { return a->get_start() < b->get_start(); });

    vector<System*> systems = gen_systems(loci);

    for (System* system : systems)
    {
        out_gene << system->to_string_summary();
         out_gene_debug << system->to_string_debug() << endl;
    }

    for (Crispr* c : crisprs) delete c;
    for (Translation* t : translations) delete t;
//    for (MultiFragment* m : multifragments) delete m;
    for (System* s : systems) delete s;

    auto timed_prospect = time_diff(start_prospect, time());

    out_gene << fmt::format("// finished in {} ms", timed_prospect);

    out_gene.close();
    out_gene_debug.close();
}

void run(vector<CasProfile*>& profiles)
{
    // interest "GCA_002139875.1_ASM213987v1_genomic.fna"
    unordered_set<string> interest {  };
    ui limiter = 1;
    ui track = 0;
    for (const auto& entry : std::filesystem::directory_iterator(Config::path_genome))
    {
        string filename = entry.path().filename().string();

        if (interest.empty() || (!interest.empty() && interest.contains(filename)))
        {
            prospect_genome(profiles, entry);
        }

        if (++track == limiter) {
            break;
        }
    }
}

struct Hit
{
    CasProfile* profile;
    ui hit_count;
};

// void hmm_db_experiment()
// {
//     vector<string> hmm_files = get_all_hmm_files();

//     std::filesystem::path path_prodigal_proteins = "/home/ben/crispr/prospector-util/my.proteins.faa";

//     map<string, string> proteins = Util::parse_fasta(path_prodigal_proteins, false);

//     map<string, vector<Hit*>> protein_id_to_hits;

//     for (auto & [identifier, sequence] : proteins)
//     {
// 		sequence.erase(remove(sequence.begin(), sequence.end(), '*'), sequence.end());
//     }

//     for (auto & [identifier, sequence] : proteins)
//     {
//         fmt::print("{}\n", identifier);
//         for (string& hmm_file : hmm_files)F
//         {
//             bool positive = singleton_hammer(sequence, hmm_file);
//             if (positive) {
//                 fmt::print("got a positive on {}\n", sequence);
//             }
//         }
//     }

// }

bool claim_validation(string& protein_sequence, string cas_claim)
{
    vector<string> hmm_files = get_hmm_files(cas_claim);
    fmt::print("executing {} hmm files\n", hmm_files.size());
    for (string hmm_file : hmm_files)
    {
        if (singleton_hammer(protein_sequence, hmm_file))
        {
            return true;
        }
    }
    return false;
}


string hit_assessment(Hit* hit, string& protein_sequence)
{
    string cas_claim = CasProfileUtil::domain_table_fetch(hit->profile->identifier);
    bool validated = claim_validation(protein_sequence, cas_claim);
    return fmt::format("{} {} {} {}\n", hit->profile->identifier, cas_claim, hit->hit_count, validated);
    
}

void analyze_prodigal_proteins(vector<CasProfile*>& profiles)
{
    auto start_prodigal = time();

    std::filesystem::path path_prodigal_proteins = "/home/ben/crispr/prospector-util/my.proteins.faa";

    map<string, string> proteins = Util::parse_fasta(path_prodigal_proteins, false);

    map<string, vector<Hit*>> protein_id_to_hits;

    for (auto & [identifier, sequence] : proteins)
    {
		sequence.erase(remove(sequence.begin(), sequence.end(), '*'), sequence.end());

        auto kmers = Util::kmerize(sequence, 6);
        auto kmers_encoded = Util::encode_amino_kmers(kmers, 6);

        vector<Hit*> hits;

        for (signed long p = 0; p < profiles.size(); p++)
        {
            CasProfile* profile = profiles[p];
            // vector<ull> index;
            ull count = 0;
            for (ull i = 0; i < kmers_encoded.size(); i++)
            {
                bool contains = profile->hash_table.contains(kmers_encoded[i]);
                if (contains)
                    count++;
                    // index.push_back(i);
            }

            Hit* hit = new Hit;
            hit->profile = profile;
            // hit->hit_count = index.size();
            hit->hit_count = count;

            hits.push_back(hit);
        }

        // sort hits
        Util::sort(hits, [](Hit* a, Hit* b) { return a->hit_count > b->hit_count; });

        // store hits
        protein_id_to_hits[identifier] = hits;
    }

    start_prodigal = time(start_prodigal, "hit counts");

    std::ofstream outtie("outtie.txt");

    for (auto const& [identifier, hits] : protein_id_to_hits)
    {
        string protein_sequence = proteins[identifier];
        outtie << fmt::format("{}\n", identifier);
        outtie << hit_assessment(hits[0], protein_sequence);
        outtie << hit_assessment(hits[1], protein_sequence);
        outtie << hit_assessment(hits[2], protein_sequence);
        outtie << "\n\n";
    }

    start_prodigal = time(start_prodigal, "hmmer");

    outtie.close();

    start_prodigal = time(start_prodigal, "full prodigal analysis");
}

int main(int argc, char *argv[])
{
    Config::parse_program_args(argc, argv);

    auto start_main = time();

    CasProfileUtil::load_domain_map(Config::path_map_dom);

    if (Config::skip_serialisation == 0) CasProfileUtil::serialize(Config::path_bin_pro);
    vector<CasProfile*> profiles = CasProfileUtil::deserialize_profiles(Config::path_bin_pro);
    Prospector::device_init();

    run(profiles);

    for (CasProfile* p : profiles) delete p;

    start_main = time(start_main, "main");
    return 0;                                                                                                           
}
