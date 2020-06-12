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


bool Array::mutant(const char* genome, ui* genome_encoding, ui k, ui i, ui j)
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

     //for (ui __i = checked_so_far; i < k; __i++)
     //{
     //    diff += genome[i + checked_so_far + __i] == genome[j + checked_so_far + __i] ? 0 : 1; 
     //}
     //return diff <= allowed_mutations;

}

vector<Crispr*> prospector_main(string& genome)
{
    Prospector::Encoding encoding = Prospector::get_genome_encoding(genome.c_str(), genome.size());
    uc* qmap = Prospector::get_qmap_small(encoding.d, encoding.size);
    vector<ui> query_indices = get_candidate_queries(qmap, encoding.size);

    map<ui, bool> consumed; for (ui query : query_indices) consumed[query] = false;

    vector<Crispr*> all_crisprs;
    for (ui query : query_indices)
    {
        if (consumed[query])
            continue;

        //if (query != 2110970)
            //continue;

        for (ui k = Prospector::k_start; k < Prospector::k_end; k++) 
        {
            vector<ui> proximals{ query };
            ui target = query + k + Prospector::spacer_min;

            while (target - ( proximals[proximals.size()-1] + k) <= Prospector::spacer_max)
            {
                target++;

                auto diff = Util::difference_cpu(encoding.h[query], encoding.h[target]);
                auto query_str = genome.substr(query, k);
                auto target_str = genome.substr(target, k);
                auto mutant = Array::mutant(genome.c_str(), encoding.h, k, query, target);


                //if (mutant)
                    //fmt::print("{} -> {} {} {} {} {}\n", query_str, target, target_str, target + k, diff, mutant);

                if (Array::mutant(genome.c_str(), encoding.h, k, query, target))
                {
                    consumed[target] = true;
                    proximals.push_back(target);
                    target += k + Prospector::spacer_min;
                }
            }

            //break;

            if (proximals.size() >= Prospector::repeats_min)
            {
                Crispr* crispr = new Crispr(k, proximals, proximals.size());
                all_crisprs.push_back(crispr);
            }
        }
    }
    return all_crisprs;
}


vector<Crispr*> Array::get_crisprs(string& genome)
{
    vector<Crispr*> crisprs = prospector_main(genome);
    CrisprUtil::cache_crispr_information(genome, crisprs);

    //Debug::crispr_print(crisprs, genome, 2111914 - 5000, 2112599 + 5000);

    crisprs = Util::filter(crisprs, [](Crispr* c) { return c->overall_heuristic >= -1.5; });
    Util::sort(crisprs, CrisprUtil::heuristic_greater);

    //Debug::crispr_print(crisprs, genome, 2643367 - 100, 2661244 + 100);

    crisprs = CrisprUtil::get_domain_best(crisprs);
    Util::sort(crisprs, [](Crispr* a, Crispr* b) { return a->start < b->start; });
    return crisprs;
}
