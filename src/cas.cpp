#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "cas.h"
#include "config.h"

static unordered_set<string> alternate_starts_pos{"GTG", "TTG"};
static unordered_set<string> alternate_starts_neg{"CAC", "CAA"};

static unordered_set<string> start_codons_pos{"ATG"};
static unordered_set<string> start_codons_neg{"CAT"};
static unordered_set<string> stop_codons_pos{"TAA", "TAG", "TGA"};
static unordered_set<string> stop_codons_neg{"TTA", "CTA", "TCA"};

// potential performance optimization: compare amino_family acids of codons in the amino_family rather than codons in the genome.
// should be faster to compare single chars rather than strings.

ull generic_expand(const string &genome, unordered_set<string> &sinks, ull begin, ull increment, ull increments) {
    ull offset = 0;
    try {
        for (ull i = 0; i < increments; i++) {
            if (sinks.contains(genome.substr(begin + offset, 3))) {
                return begin + offset;
            }
            offset += increment;
        }
    }
    catch (exception &) {
        return begin;
    }
    return begin;
}

ull steps = 200;

void fragment_expansion_pos(Fragment *fragment, const string &genome) {
    fragment->expanded_genome_begin = generic_expand(genome, start_codons_pos, fragment->genome_begin, -3, steps);
    fragment->expanded_genome_final = generic_expand(genome, stop_codons_pos, fragment->genome_final - 3, 3, steps) + 3;

    if (fragment->expanded_genome_begin == fragment->genome_begin) {
        fragment->expanded_genome_begin = generic_expand(genome, alternate_starts_pos, fragment->genome_begin, -3,
                                                         steps);
    }
}

void fragment_expansion_neg(Fragment *fragment, const string &genome) {
    fragment->expanded_genome_begin = generic_expand(genome, stop_codons_neg, fragment->genome_begin, -3, steps);
    fragment->expanded_genome_final =
            generic_expand(genome, start_codons_neg, fragment->genome_final - 3, 3, steps) + 3;

    if (fragment->expanded_genome_final == fragment->genome_final) {
        fragment->expanded_genome_final =
                generic_expand(genome, alternate_starts_neg, fragment->genome_final - 3, 3, steps) + 3;
    }

}

vector<vector<ull>> cluster_index(const vector<ull> &indices) {
    vector<vector<ull>> clusters;
    vector<ull> cluster;
    ull prev = indices[0];
    for (ull index: indices) {
        if (index - prev > Cas::max_inter_cluster_dist) {
            vector<ull> cluster_cp = cluster;
            if (cluster_cp.size() > 1) {
                clusters.push_back(cluster_cp);
            }
            cluster.clear();
        }
        cluster.push_back(index);
        prev = index;
    }
    clusters.push_back(cluster);
    return clusters;
}

bool good_clusters(const vector<vector<ull>> &clusters) {
    return std::ranges::any_of(clusters.begin(), clusters.end(),
                               [](auto &cluster) { return cluster.size() > Cas::cluster_metric_min; });
}

ull demarc_start_clusters(const vector<vector<ull>> &clusters) {
    for (const vector<ull> &cluster: clusters)
        if (cluster.size() > Cas::cluster_metric_min)
            return cluster[0];
    assert(false);
    return -1;
}

ull demarc_final_clusters(const vector<vector<ull>> &clusters) {
    for (ull i = clusters.size() - 1; i >= 0; i--)
        if (clusters[i].size() > Cas::cluster_metric_min)
            return clusters[i][clusters[i].size() - 1];
    assert(false);
    return -1;
}

void compute_demarc(Fragment *frag) {
    frag->clust_begin = demarc_start_clusters(frag->clusters);
    frag->clust_final = demarc_final_clusters(frag->clusters);
}

void compute_details(Fragment *fragment, const string &genome) {
    ull index_kmer_start = fragment->clust_begin;
    ull index_kmer_final = fragment->clust_final;

    fragment->protein_sequence = fragment->reference_translation->pure.substr(index_kmer_start,
                                                                              (index_kmer_final - index_kmer_start) +
                                                                              CasProfileUtil::k);

    ull raw_pos_start = fragment->reference_translation->pure_mapping[index_kmer_start];
    ull raw_pos_end = fragment->reference_translation->pure_mapping[index_kmer_final];

    ull genome_begin;
    ull genome_final;
    if (fragment->reference_translation->pos) {
        genome_begin = fragment->reference_translation->genome_start + (raw_pos_start * 3); // inclusive begin
        genome_final = fragment->reference_translation->genome_start +
                       ((raw_pos_end + CasProfileUtil::k) * 3); // exclusive final
    } else {
        genome_final = fragment->reference_translation->genome_final - (raw_pos_start * 3); // exclusive final
        genome_begin = fragment->reference_translation->genome_final -
                       (((raw_pos_end + CasProfileUtil::k) * 3)); // inclusive begin
    }

    genome_final += 1; // to be exclusive.

    //string genome_translation = Debug::translation_test(genome, genome_start, genome_final, f.reference_translation->pos, 6);
    //string amino_family = Util::translate_genome(genome, genome_start, genome_final, f.reference_translation->pos);


    fragment->genome_begin = genome_begin;
    fragment->genome_final = genome_final;
    //f->genome_translation = genome_translation;
    //details->quality = genome_final - genome_start;

    //f.details = details;
}

bool fragment_equivalent(const Fragment *a, const Fragment *b) {
    return a->clust_begin == b->clust_begin && a->clust_final == b->clust_final;
}

bool fragment_contains(const Fragment *a, const Fragment *b) {
    return a->clust_begin > b->clust_begin && a->clust_final < b->clust_final;
}

vector<Prediction *> collapse_predictions(vector<Prediction *> predictions) {
    vector<Prediction *> collapsed_predictions;
    if (predictions.empty()) return collapsed_predictions;

    // Sort predictions by start bp, ascending
    std::sort(predictions.begin(), predictions.end(), [](auto &a, auto &b) { return a->start < b->start; });

    int skip;
    for (int i = 0; i < predictions.size() - 1; i += skip) {
        skip = 1;
        auto p1 = predictions[i];
        auto p2 = predictions[i + 1];
        double total_score = p1->score;

        while (p1->id == p2->id && (i + skip + 1) < predictions.size()) {
            // Both predictions have the same Cas family (but potentially different profiles)
            // They can be 'collapsed' into a single prediction over a given domain
            p1 = predictions[i + skip];
            p2 = predictions[i + skip + 1];

            total_score += p1->score;
            skip++;
        }

        auto *p = new Prediction;
        p->score = total_score / skip;
        p->start = predictions[i]->start;
        p->end = predictions[i + skip - 1]->end;
        p->id = predictions[i]->id;

        collapsed_predictions.push_back(p);
    }

    return collapsed_predictions;
}

vector<Prediction *> Cas::predict_cas(vector<CasProfile *> &profiles, string &genome, ull offset) {
    // Load frequency map
    FrequencyMap frequency_map;
    frequency_map.load();

    phmap::flat_hash_map<size_t, Prediction *> predicted_regions_map;
    size_t iterations = genome.size() / Config::cas_chunk_length;

#pragma omp parallel for
    for (size_t i = 0; i < iterations; i++) {
        size_t start = i * Config::cas_chunk_length;
        size_t end = start + Config::cas_chunk_length;

        // Sample a chunk of the genome at a time
        string genome_substr = genome.substr(start, Config::cas_chunk_length);
        vector<Translation *> translations = Cas::get_sixframe(genome_substr, 0, genome_substr.length() - 1);

        // For each profile, count the number of (weighted) kmers in the genome chunk
        for (auto &profile: profiles) {
            double score = 0.0;

            for (auto &t: translations)
                for (auto &km: t->pure_kmerized_encoded)
                    if (profile->hash_table.contains(km)) score += frequency_map.freq_map[km];

#pragma omp critical
            {
                bool should_save_score =
                        !predicted_regions_map.contains(start) || predicted_regions_map[start]->score < score;

                // Store predictions by the region they started in
                if (score > Config::cas_threshold && should_save_score) {
                    auto *p = new Prediction;
                    p->score = score;
                    p->start = start + offset;
                    p->end = end + offset;
                    p->id = CasProfileUtil::domain_table_fetch(profile->identifier);

                    predicted_regions_map[start] = p;
                }
            };
        }
    }

    // Convert from map to a list
    vector<Prediction *> predictions;
    for (auto &p: predicted_regions_map) predictions.push_back(p.second);

    // Merge chunk predictions into larger chunks
    auto collapsed_pred = collapse_predictions(predictions);

    // Cleanup old predictions
    for (Prediction* p : predictions) delete p;

    return collapsed_pred;
}

vector<Fragment *> old_cas(vector<CasProfile *> &profiles, vector<Translation *> &translations, string &genome) {
    vector<Fragment *> fragments;

#pragma omp parallel for
    for (auto profile: profiles) {
        ui translation_index = -1;
        for (Translation *t: translations) {
            translation_index++;
            vector<ull> index;

            for (ull i = 0; i < t->pure_kmerized_encoded.size(); i++) {
                bool contains = profile->hash_table.contains(t->pure_kmerized_encoded[i]);
                if (contains)
                    index.push_back(i);
            }

            if (index.size() == 0)
                continue;

            // dump index to a file for analysis
            string file_name = Config::path_output / fmt::format("index_dump/{}_{}_{}_{}_{}_{}",
                                                                 CasProfileUtil::domain_table_fetch(
                                                                         profile->identifier), profile->identifier,
                                                                 t->genome_start, t->genome_final, translation_index,
                                                                 index.size());

            if (Config::dump_indices) {
                std::ofstream out(file_name);
                for (ull index_location: index) {
                    out << index_location << "\t" << t->genome_start + index_location << endl;
                }
                out.close();
            }

            vector<vector<ull>> clusters = cluster_index(index);
            if (!good_clusters(clusters)) {
                continue;
            } else if (Config::dump_indices) {
                fmt::print("accepted {}\n", file_name);

                string file_name_clust = fmt::format("{}.clust", file_name);
                std::ofstream out(file_name_clust);
                for (const vector<ull> &clust: clusters) {
                    for (ull index_location: clust) {
                        out << index_location << "\t" << t->genome_start + index_location << endl;
                    }
                    out << endl;
                }
                out.close();
            }


            auto *f = new Fragment;
            f->reference_genome = &genome;
            f->reference_crispr = t->reference_crispr;
            f->reference_translation = t;
            f->reference_profile = profile;
            f->clusters = clusters;

            compute_demarc(f);

            if (f->clust_final - f->clust_begin <= 15)
                continue;

            compute_details(f, genome);

            auto expansion = f->reference_translation->pos ? fragment_expansion_pos : fragment_expansion_neg;
            expansion(f, genome);

            if (f->genome_begin < 0) {
                std::exit(1);
            }

            // if (f->repetition_problem())
            // {
            // fmt::print("rejected on the basis of repetition problem\n");
            // continue;
            // }

#pragma omp critical
            {
                fragments.push_back(f);
            }
        }
    }

    std::sort(fragments.begin(), fragments.end(),
              [](Fragment *a, Fragment *b) { return a->expanded_genome_begin < b->expanded_genome_begin; });

    return fragments;
}

vector<Translation *> Cas::get_triframe(const string &genome, ull genome_start, ull genome_final, bool pos) {
    string domain = genome.substr(genome_start, genome_final - genome_start);
    //domain = pos ? domain : Util::reverse_complement(domain);

    if (!pos) {
        Util::reverse_complement(domain);
    }

    vector<Translation *> translations;
    for (ull frame = 0; frame < 3; frame++) {
        auto *translation = new Translation;
        translation->pos = pos;
        translation->genome_start = pos ? genome_start + frame : genome_start;
        translation->genome_final = pos ? genome_final : genome_final - frame;

        translation->raw = Util::translate_domain(domain.substr(frame));

        translation->pure = "";

        ull stop_count = 0;
        ull index = 0;
        for (char elem: translation->raw) {
            if (elem == Util::stop_c) {
                stop_count++;
                continue;
            }
            translation->pure += elem;
            translation->pure_mapping.push_back(index + stop_count);
            index++;
        }
        translation->pure_kmerized = Util::kmerize(translation->pure, CasProfileUtil::k);
        translation->pure_kmerized_encoded = Util::encode_amino_kmers(translation->pure_kmerized, CasProfileUtil::k);
        translations.push_back(translation);
    }
    return translations;
}

vector<Translation *> Cas::get_sixframe(const string &genome, ull genome_start, ull genome_final) {
    vector<Translation *> sixframe;

    for (Translation *translation: Cas::get_triframe(genome, genome_start, genome_final, true))
        sixframe.push_back(translation);

    for (Translation *translation: Cas::get_triframe(genome, genome_start, genome_final, false))
        sixframe.push_back(translation);

    return sixframe;
}

vector<Translation *> Cas::crispr_proximal_translations(const string &genome, vector<Crispr *> &crisprs) {
    vector<Translation *> translations;

    std::sort(crisprs.begin(), crisprs.end(), [](Crispr *a, Crispr *b) { return a->genome_start < b->genome_start; });
    const ull min_translation_size = 20;
    const ull permitted_translation_overlap = 100;
    ull minimum_bound = 0;
    for (Crispr *c: crisprs) {
        ull translation_start;
        ull translation_final;

        // normal case
        translation_start = c->genome_start - Cas::upstream_size;
        translation_final = c->genome_final + Cas::upstream_size;

        // guard bounds
        translation_start = translation_start > c->genome_start ? 0 : translation_start; // underflow
        translation_final = std::min(translation_final, (ull) genome.size() - 1); // exceeds genome

        // bound the start
        translation_start = std::max(translation_start, minimum_bound);

        // require min size
        if (translation_final - translation_start < min_translation_size)
            continue;

        // new bound
        minimum_bound = translation_final - permitted_translation_overlap;
        assert(minimum_bound < translation_final); // guard underflow

        vector<Translation *> six_frame = Cas::get_sixframe(genome, translation_start, translation_final);
        for (Translation *t: six_frame) t->reference_crispr = c;
        translations.insert(translations.end(), six_frame.begin(), six_frame.end());
    }
    return translations;
}


#pragma clang diagnostic pop