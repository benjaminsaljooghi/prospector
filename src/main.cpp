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
#include "system.h"

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
    auto start_prospect = time();

    fmt::print("\n\n");

    std::filesystem::path results_path = Config::path_results / genome_path.stem();    
    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_gene_debug(results_path / "out_gene_debug.txt");

    string genome = Util::load_genome(genome_path);

    vector<Translation*> translations;
    vector<Crispr*> crisprs;

    if (Config::cas_only) {
        translations = Cas::get_sixframe(genome, 0, genome.length()-1);
    } else {
        crisprs = Array::get_crisprs(genome);
        translations = Config::crispr_proximal_search ? Cas::crispr_proximal_translations(genome, crisprs) : Cas::get_sixframe(genome, 0, genome.length()-1);
    }

    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
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
        out_gene_debug << system->to_string_debug() << endl;
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

void run(vector<CasProfile*>& profiles)
{
    unordered_set<string> interest {};
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
