#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time_util.h"
#include "cas_profiles.h"
#include "array_discovery.h"
#include "path.h"
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

    if (loci.empty())
        return systems;

    auto* current = new System;
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

    Util::GenomeIdSequenceMap genome_map = Util::load_genome(genome_path);

    for (auto const& [genome_id, genome_sequence] : genome_map)
    {
        vector<Crispr*> crisprs;
        vector<Prediction*> predictions;

        fmt::print("Acquiring CRISPRs & translations for {}...\n", genome_id);
        crisprs = Array::get_crisprs(const_cast<string &>(genome_sequence));

        ull min_start = ULONG_LONG_MAX; ull max_end = 0;
        for (auto& c: crisprs) {
            min_start = min(min_start, c->genome_start);
            max_end = max(max_end, c->genome_final);
        }

        string cas_region = genome_sequence.substr(min_start - Cas::upstream_size, (max_end - min_start) + Cas::upstream_size);

        fmt::print("Detecting Cas genes for {}...\n", genome_id);
        predictions = Cas::predict_cas(profiles, const_cast<string &>(cas_region), min_start - Cas::upstream_size);

        fmt::print("Collating results for {}...\n", genome_id);
        std::vector<Locus*> loci;

        for (Crispr* c : crisprs)
            loci.push_back(c);

        for (Prediction* p: predictions)
            loci.push_back(p);

        std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) { return a->get_start() < b->get_start(); });

        vector<System*> systems = gen_systems(loci);

        fmt::print("Writing results for {} to file...\n", genome_id);
        for (System* system : systems)
        {
            out_gene << system->to_string_summary(const_cast<string &>(genome_id));
            out_gene_debug << system->to_string_debug(const_cast<string &>(genome_id)) << endl;
        }

        for (Crispr* c : crisprs) delete c;
        for (Prediction* p : predictions) delete p;
        for (System* s : systems) delete s;
    }

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
