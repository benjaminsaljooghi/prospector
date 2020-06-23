#include "array_discovery.h"


vector<ui> get_candidate_queries(unsigned char* qmap, ui genome_encoding_size)
{
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
    return query_indices;
}


bool Array::mutant(const char* genome, ui* genome_encoding, ui k, ui i, ui j, ui repeat_tolerance_ratio)
{
    ui allowed_mutations = k / repeat_tolerance_ratio;

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

    return diff <= checked_so_far / repeat_tolerance_ratio;

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

        for (ui k = Prospector::k_start; k < Prospector::k_end; k++) 
        {
            vector<ui> proximals{ query };
            ui target = query + k + Prospector::spacer_min;

            while (++target - ( proximals[proximals.size()-1] + k) <= Prospector::spacer_max && target + k < genome.size())
            {
                
                //auto diff = Util::difference_cpu(encoding.h[query], encoding.h[target]);
                //auto query_str = genome.substr(query, k);
                //auto target_str = genome.substr(target, k);
                //auto mutant = Array::mutant(genome.c_str(), encoding.h, k, query, target);
                //if (mutant)
                    //fmt::print("{} -> {} {} {} {} {}\n", query_str, target, target_str, target + k, diff, mutant);

                if (Array::mutant(genome.c_str(), encoding.h, k, query, target, Prospector::repeat_tolerance_ratio))
                {
                    consumed[target] = true;
                    proximals.push_back(target);
                    target += k + Prospector::spacer_min;
                }
            }

            if (proximals.size() >= Prospector::repeats_min && proximals.size() < 500 && proximals[proximals.size()-1] > proximals[0])
            {
                // perform expansion
                ui __query;
                ui __target;

                // forward
                __query = proximals[proximals.size() - 1];
                __target = __query + k + Prospector::spacer_min;
                while (++__target - (proximals[proximals.size() - 1] + k) <= Prospector::spacer_max && target + k < genome.size())
                {
                    if (Array::mutant(genome.c_str(), encoding.h, k, __query, __target, Prospector::repeat_tolerance_ratio_sensitive))
                    {
                        proximals.push_back(__target);
                        __target += k + Prospector::spacer_min;
                    }
                }

                //backward
                __query = proximals[0];
                __target = query - Prospector::spacer_min - k;
                while (proximals[0] - (--__target + k) <= Prospector::spacer_max && __target < __query) // guard overflow
                {
                    if (Array::mutant(genome.c_str(), encoding.h, k, __query, __target, Prospector::repeat_tolerance_ratio_sensitive))
                    {
                        proximals.insert(proximals.begin(), __target);
                        __target -= (k + Prospector::spacer_min);
                    }
                }
                
                Crispr* crispr = new Crispr(k, proximals, proximals.size());
                
               
                all_crisprs.push_back(crispr);
            }
        }
    }

    Prospector::free_encoding(encoding);

    return all_crisprs;
}


vector<Crispr*> Array::get_crisprs(string& genome)
{
    auto start = time();

    vector<Crispr*> crisprs = prospector_main(genome);
    start = time(start, "prospect");
    
    CrisprUtil::cache_crispr_information(genome, crisprs);
    start = time(start, "cache");

    crisprs = Util::filter(crisprs, [](Crispr* c) { return c->overall_heuristic >= -3; });
    Util::sort(crisprs, CrisprUtil::heuristic_greater);
    crisprs = CrisprUtil::get_domain_best(crisprs);
    Util::sort(crisprs, [](Crispr* a, Crispr* b) { return a->start < b->start; });
    start = time(start, "remain");

    return crisprs;
}
