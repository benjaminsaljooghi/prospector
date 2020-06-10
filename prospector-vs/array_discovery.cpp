#include "array_discovery.h"


vector<ui> get_candidate_queries(unsigned char* qmap, ui genome_encoding_size)
{
    auto start = time();
    vector<ui> query_indices;
    static const ui tolerance = Prospector::size / Prospector::repeat_tolerance_ratio;
    for (ui query_index = 0; query_index < genome_encoding_size - 200; query_index++)
    {
        ui ting = query_index * Prospector::map_size_small;
        for (ui i = 0; i < Prospector::map_size_small; i++)
        {
            if (qmap[ting + i] <= tolerance)
            {
                query_indices.push_back(query_index);
                break;
            }
        }
    }
    time(start, "post qmap small candidate query genertion");
    return query_indices;
}


bool mutant(const char* genome, const ui* genome_encoding, const ui& k, const ui& i, const ui& j)
{
    ui allowed_mutations = k / Prospector::repeat_tolerance_ratio;

    ui diff = 0;
    const ui chunks = k / Prospector::size;
    // may generate a lot of crisprs that are filtered later (expensive) given that SIZE is large (16) here.
    // option is to do a more accurate mutation calculation either using per-char post the chunk division
    // or to encode entire kmers up to MAP_SIZ into ull's to compute the difference efficiently.
    // post k=MAP_SIZ we can use the MAP_SIZ-long ull initially, and then compute a per-char difference afterwards.

    for (ui chunk = 0; chunk < chunks; chunk++)
    {
        ui _i = genome_encoding[i + (chunk * Prospector::size)];
        ui _j = genome_encoding[j + (chunk * Prospector::size)];
        diff += Util::difference_cpu(_i, _j);
        if (diff > allowed_mutations)
        {
            return false;
        }
    }
    const ui checked_so_far = (chunks * Prospector::size);

    return diff <= checked_so_far / Prospector::repeat_tolerance_ratio;

    // for (ui __i = checked_so_far; i < k; __i++)
    // {
        // diff += genome[i + checked_so_far + __i] == genome[j + checked_so_far + __i] ? 0 : 1; 
    // }
    // return diff <= allowed_mutations;

}

vector<Crispr*> prospector_main(string& genome)
{
    Prospector::Encoding encoding = Prospector::get_genome_encoding(genome.c_str(), genome.size());
    uc* qmap = Prospector::get_qmap_small(encoding.encoding_d, encoding.size);

    //Debug::visualize_map(encoding, genome);

    vector<ui> query_indices = get_candidate_queries(qmap, encoding.size);

    map<ui, bool> consumed;

    for (ui query_index : query_indices)
        consumed[query_index] = false;

    vector<Crispr*> all_crisprs;

    // can optimize this by performing a spacer_min leap
    // can optimize this by starting at largest k, then breaking this loop on the first time we fail to form a crispr

    //ui advancement = 0
    for (ui query_index : query_indices)
    {
        if (consumed[query_index])
        {
            continue;
        }

        for (ui k = Prospector::k_start; k < Prospector::k_end; k++) 
        {
            vector<ui> proximals{ query_index };

            ui target_index = query_index + k + Prospector::spacer_min;

            while (target_index - proximals[proximals.size()-1] <= Prospector::spacer_max)
            {
                //if (query_index == 2643432 && target_index == 2657405 && k == 30)
                //{
                //    printf("break\n");
                //}

                target_index++;
                if (mutant(genome.c_str(), encoding.encoding, k, query_index, target_index))
                {
                    consumed[target_index] = true;
                    proximals.push_back(target_index);
                    target_index += k + Prospector::spacer_min;
                }
            }

            if (proximals.size() >= Prospector::repeats_min)
            {
                Crispr* crispr = new Crispr(k, proximals, proximals.size());
                all_crisprs.push_back(crispr);
            }
        }


    }

    //uc* qmap_big = Prospector::get_qmap_big(encoding.encoding_d, encoding.size, &query_indices[0], query_indices.size());

    // vector<vector<ui>> proximal_targets;
    //auto start = time();
    //map<ui, vector<ui>> proximal_targets;
    //ui tolerance = 16 / Prospector::repeat_tolerance_ratio;

    //for (ui q_i = 0; q_i < query_indices.size(); q_i++)
    //{
    //    ui query_index = query_indices[q_i];
    //    vector<ui> proximals;
    //    for (ui t_i = 0; t_i < Prospector::map_size_big; t_i++)
    //    {
    //        ui qmap_index = q_i * Prospector::map_size_big + t_i;
    //        if (qmap_big[qmap_index] <= tolerance)
    //        {
    //            ui target = query_index + Prospector::k_start + Prospector::spacer_skip + t_i;
    //            proximals.push_back(target);
    //        }
    //    }

    //    // remove any exact subsetting
    //    auto check = [&]() {
    //        for (auto const& [query_index, proximal] : proximal_targets)
    //            for (ui p : proximals)
    //                if (Util::contains(proximal, p))
    //                    return false;
    //        return true;
    //    };

    //    //if (check())
    //        proximal_targets[query_index] = proximals;

    //}


    //Debug::visualize_proximals(proximal_targets, Prospector::repeats_min, genome);

    //vector<Crispr*> all_crisprs;
    //for (auto const& [query_index, proximal] : proximal_targets)
    //{
    //    vector<Crispr*> candidates;
    //    for (ui k = Prospector::k_end - 1; k >= Prospector::k_start; k--)
    //    {
    //        ui allowed_mutations = k / Prospector::repeat_tolerance_ratio;
    //        vector<ui> genome_indices;
    //        genome_indices.push_back(query_index);
    //        for (ui target : proximal)
    //        {
    //            ui end = genome_indices[genome_indices.size() - 1] + k;
    //            if (target < end || target - end < Prospector::spacer_min) continue; // || guards against overflow
    //            if (target - end > Prospector::spacer_max) break;

    //            if (mutant(genome.c_str(), encoding.encoding, k, allowed_mutations, query_index, target))
    //                genome_indices.push_back(target);
    //        }

    //        if (genome_indices.size() >= Prospector::repeats_min)
    //        {
    //            Crispr* c = new Crispr(k, genome_indices, genome_indices.size());
    //            candidates.push_back(c);
    //        }
    //    }

    //    if (candidates.size() > 0)
    //    {
    //        Util::sort(candidates, [](const Crispr* a, const Crispr* b) {  return a->size > b->size; });
    //        ui best_size = candidates[0]->size;
    //        for (Crispr* crispr : candidates)
    //        {
    //            if (crispr->size != best_size)
    //                break;

    //            all_crisprs.push_back(crispr);
    //        }
    //    }
    //}




    //time(start, "post-kernel crispr generation");
    //fmt::print("\tprospector returned {} crisprs\n", all_crisprs.size());
    //return all_crisprs;
    return all_crisprs;
}


vector<Crispr*> Array::get_crisprs(string& genome)
{
    vector<Crispr*> crisprs = prospector_main(genome);
    CrisprUtil::cache_crispr_information(genome, crisprs);

    crisprs = Util::filter(crisprs, [](Crispr* c) { return c->overall_heuristic >= -2; });
    Util::sort(crisprs, CrisprUtil::heuristic_greater);

    //Debug::crispr_print(crisprs, genome, 2643367 - 100, 2661244 + 100);

    crisprs = CrisprUtil::get_domain_best(crisprs);
    Util::sort(crisprs, [](Crispr* a, Crispr* b) { return a->start < b->start; });
    return crisprs;
}
