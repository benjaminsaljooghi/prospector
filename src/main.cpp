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

#define MAX_MISMATCHES 500


vector<MultiFragment *> gen_multifragments(vector<Fragment *> fragments) {
    vector < MultiFragment * > multifragments;
    for (ui i = 0; i < fragments.size(); i++) {
        fmt::print("{}: {}\n", i, CasProfileUtil::domain_table_fetch(fragments[i]->reference_profile->identifier));

        fmt::print("multifragment {}\n", i);
        MultiFragment *multifragment = new MultiFragment;
        multifragment->fragments.push_back(fragments[i]);

        for (ui j = i + 1; j < fragments.size(); j++) {
            bool any_overlap = Util::any_overlap(fragments[i]->genome_begin, fragments[i]->genome_final,
                                                 fragments[j]->genome_begin, fragments[j]->genome_final);

            string first = CasProfileUtil::domain_table_fetch(fragments[i]->reference_profile->identifier);
            string second = CasProfileUtil::domain_table_fetch(fragments[j]->reference_profile->identifier);

            bool domain_overlap = (first.find(second) != string::npos) || (second.find(first) != string::npos);

            if (any_overlap && domain_overlap) {
                multifragment->fragments.push_back(fragments[j]);
                i = j;
            }
        }

        multifragments.push_back(multifragment);
    }
    return multifragments;
}

vector<System *> gen_systems(vector<Locus *> loci) {
    vector < System * > systems;

    if (loci.size() == 0)
        return systems;

    System *current = new System;
    current->loci.push_back(loci[0]);
    for (size_t i = 1; i < loci.size(); i++) {
        if (loci[i]->get_start() < current->get_final() + 30000) {
            current->loci.push_back(loci[i]);
        } else {
            systems.push_back(current);
            current = new System;
            current->loci.push_back(loci[i]);
        }
    }

    systems.push_back(current);

    vector < System * > filtered_systems;

    for (System *s: systems) {
        if (s->legitimate_system()) {
            filtered_systems.push_back(s);
        } else {
            delete s;
        }
    }

    return filtered_systems;
}

void print_int_kmer(uint64_t k) {
    for (uint64_t i = 0; i < 8 * CasProfileUtil::k; i += 8) {
        uint64_t offset = (uint64_t) 0xFF << i;
        uint64_t c = ((k & offset) >> i);
        printf("%c", (char) c);
    }

    printf("\n");
}

void prospect_genome(vector<CasProfile *> &profiles, std::filesystem::path genome_path) {
    auto start_prospect = time();

    fmt::print("\n\n");

    std::filesystem::path results_path = Config::path_results / genome_path.stem();
    std::filesystem::create_directory(results_path);
    std::ofstream out_gene(results_path / "out_gene.txt");
    std::ofstream out_gene_debug(results_path / "out_gene_debug.txt");

    Util::GenomeIdSequenceMap genome_map = Util::load_genome(genome_path);

    for (auto const &[genome_id, genome_sequence]: genome_map) {
        vector < Translation * > translations;
        vector < Crispr * > crisprs;

        if (Config::cas_only) {
            fmt::print("Acquiring translations for {}...\n", genome_id);
            translations = Cas::get_sixframe(genome_sequence, 0, genome_sequence.length() - 1);
        } else {
            fmt::print("Acquiring CRISPRs & translations for {}...\n", genome_id);
            crisprs = Array::get_crisprs(const_cast<string &>(genome_sequence));
            translations = Config::crispr_proximal_search ?
                           Cas::crispr_proximal_translations(genome_sequence, crisprs) :
                           Cas::get_sixframe(genome_sequence, 0, genome_sequence.length() - 1);
        }

        fmt::print("Detecting Cas genes for {}...\n", genome_id);

        double highest_accuracy = 0.0;
        string best_profile = "NONE";

#pragma omp parallel for
        for (auto &p: profiles) {
//            if (p->identifier == "COG3513") {
            for (auto &t: translations) {
                ui profile_size = p->binary_kmers.size();
                ui trans_size = t->binary_kmers.size();

                ui iterations = trans_size - profile_size;

                for (ui iteration = 0; iteration < iterations; iteration++) {
                    ui n_mismatch = 0;

                    for (ui i = 0; i < profile_size; i++) {
                        if (i + iteration > trans_size) break;

                        uint64_t p_k = p->binary_kmers[i];
                        uint64_t p_mask = p->binary_masks[i];
                        uint64_t t_k = t->binary_kmers[i + iteration];

                        uint64_t p_k_masked = p_k & p_mask;
                        uint64_t t_k_masked = t_k & p_mask;


                        if (p_k_masked != t_k_masked) {
                            n_mismatch++;
                        }


//                            if (n_mismatch > MAX_MISMATCHES) break;
                    }

                    double accuracy = ((double) profile_size - (double) n_mismatch) / (double) profile_size;

#pragma omp critical
                    if (accuracy > highest_accuracy) {
                        highest_accuracy = accuracy;
                        best_profile = p->identifier;
                    }

//                    if (n_mismatch < MAX_MISMATCHES) {
//                        printf("Between %u and %u, there were %u mismatches.\n",
//                               iteration,
//                               profile_size + iteration, n_mismatch);
//                    }
                }

//                fmt::print("Profile size is (in kmers): {}. Least mismatches for translation was {}...\n",
//                           profile_size, least_mismatches);
            }
//            }
        }

        fmt::print("{} had the best performance with {}% accuracy!\n", best_profile, highest_accuracy);

        break;
//        vector<Fragment*> fragments = Cas::cas(profiles, translations, const_cast<string &>(genome_sequence));
        vector < Fragment * > fragments = Cas::cas_gpu(profiles, translations, const_cast<string &>(genome_sequence));

        vector < MultiFragment * > multifragments = gen_multifragments(fragments);

        fmt::print("Collating results for {}...\n", genome_id);
        std::vector<Locus *> loci;

        for (Crispr *c: crisprs)
            loci.push_back(c);

        for (MultiFragment *f: multifragments)
            loci.push_back(f);

        std::sort(loci.begin(), loci.end(), [](Locus *a, Locus *b) { return a->get_start() < b->get_start(); });

        vector < System * > systems = gen_systems(loci);

        fmt::print("Writing results for {} to file...\n", genome_id);
        for (System *system: systems) {
            out_gene << system->to_string_summary(const_cast<string &>(genome_id));
            out_gene_debug << system->to_string_debug(const_cast<string &>(genome_id)) << endl;
        }

        for (Crispr *c: crisprs) delete c;
        for (Translation *t: translations) delete t;
        for (MultiFragment *m: multifragments) delete m;
        for (System *s: systems) delete s;
    }

    auto timed_prospect = time_diff(start_prospect, time());

    out_gene << fmt::format("// finished in {} ms", timed_prospect);

    out_gene.close();
    out_gene_debug.close();
}

void run(vector<CasProfile *> &profiles) {
    unordered_set <string> interest{};
    ui limiter = 1;
    ui track = 0;
    for (const auto &entry: std::filesystem::directory_iterator(Config::path_genome)) {
        string filename = entry.path().filename().string();

        if (interest.empty() || (!interest.empty() && interest.contains(filename))) {
            prospect_genome(profiles, entry);
        }

        if (++track == limiter) {
            break;
        }
    }
}

int main(int argc, char *argv[]) {
    Config::parse_program_args(argc, argv);

    auto start_main = time();

    CasProfileUtil::load_domain_map(Config::path_map_dom);

    if (Config::skip_serialisation == 0) CasProfileUtil::serialize(Config::path_bin_pro);
    vector < CasProfile * > profiles = CasProfileUtil::deserialize_profiles(Config::path_bin_pro);
    Prospector::device_init();

    run(profiles);

    for (CasProfile *p: profiles) delete p;

    start_main = time(start_main, "main");
    return 0;
}
