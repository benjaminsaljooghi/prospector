//std
#include "stdafx.h"

//proj
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"
#include "cas.h"

#include "time.h"


ui difference_cpu(const ui& _a, const ui& _b)
{
    ui _xor = (_a ^ _b);
    ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ui oddBits = _xor & 0x5555555555555555ull;
    ui comp = (evenBits >> 1) | oddBits;
    return __builtin_popcount(comp);
}


vector<ui> get_candidate_queries(unsigned char* qmap, ui genome_encoding_size)
{
    // how many qs in this map are containment oriented
    auto start = time();
    // ui count = 0;
    vector<ui> queries;
    for (ui query = 0; query < genome_encoding_size - 200; query++)
    {
        for (ui i = 0; i < Prospector::map_size_small; i++)
        {
            if (qmap[(query * Prospector::map_size_small) + i] <= (Prospector::size / Prospector::mutant_tolerance_ratio))
            {
                queries.push_back(query);
                break;
            }
        }
    }
    time(start, "post qmap small candidate query genertion");
    // printf("%d %zd\n", genome_encoding_size-200, queries.size());
    // return count;
    return queries;
}


bool mutant(const char* genome, const ui* genome_encoding, const ui& k, const ui& allowed_mutations, const ui& i, const ui& j)
{
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
        diff += difference_cpu(_i, _j);
        if (diff > allowed_mutations)
        {
            return false;
        }
    }
    const ui checked_so_far = (chunks * Prospector::size);

    return diff <= checked_so_far / Prospector::mutant_tolerance_ratio;

    // for (ui __i = checked_so_far; i < k; __i++)
    // {
        // diff += genome[i + checked_so_far + __i] == genome[j + checked_so_far + __i] ? 0 : 1; 
    // }
    // return diff <= allowed_mutations;
    
}


// void debug_map()
// {
//     // ui query = 1283501;
//     // ui q = genome_encoding[query];
//     // for (ui i = 0; i < 1000; i++)
//     // {
//     //     ui pos = query + Prospector::k_start + Prospector::spacer_skip + i;
//     //     ui diff = difference_cpu(genome_encoding[query], genome_encoding[pos]);

//     //     printf("%s %d %d\n", genome.substr(pos, SIZE).c_str(), pos, diff);
//     // }
// }


// visualzie proximals
// for ( auto const& [query, proximal] : proximal_targets)
// {
//     // vector<ui> proximals = proximal_targets[q_i];
//     if (proximal.size() >= Prospector::repeats_min)
//     {
//         // ui query = queries[q_i];
//         fmt::print("{}:{}-{}\n", query, genome.substr(query, 16), genome.substr(query+16, 16));
//         for (ui target : proximal)
//         {
//             // ui qmap_index = q_i * Prospector::map_size_big + t_i;
//             // ui target = query + Prospector::k_start + Prospector::spacer_skip + t_i;
//             fmt::print("\t{}:{}-{}\n", target, genome.substr(target, 16), genome.substr(target+16, 16) );
//         }
//     }
//     fmt::print("\n");
// }


vector<Crispr> prospector_main(const string& genome)
{
    Prospector::Encoding encoding = Prospector::get_genome_encoding(genome.c_str(), genome.size());

    uc* qmap = Prospector::get_qmap_small(encoding.encoding_d, encoding.size);

    vector<ui> queries = get_candidate_queries(qmap, encoding.size);
    
    uc* qmap_big = Prospector::get_qmap_big(encoding.encoding_d, encoding.size, &queries[0], queries.size());

    // vector<vector<ui>> proximal_targets;

    auto start = time();

    map<ui, vector<ui>> proximal_targets;
    ui tolerance = 16 / Prospector::mutant_tolerance_ratio;

    for (ui q_i = 0; q_i < queries.size(); q_i++)
    {
        ui query = queries[q_i];
        vector<ui> proximals;
        for (ui t_i = 0; t_i < Prospector::map_size_big; t_i++)
        {
            ui qmap_index = q_i*Prospector::map_size_big+t_i;
            if (qmap_big[qmap_index] <= tolerance)
            {
                ui target = query + Prospector::k_start + Prospector::spacer_skip + t_i;
                proximals.push_back(target);
            }
        }

        // if any of these proximals are included in an existing proximal, then we cannot include
        bool include = true;

        // remove any exact subsetting

        for ( auto const& [query, proximal] : proximal_targets)
        {
            for (ui p : proximals)
            {
                if (Util::contains(proximal, p))
                {
                    include = false;
                    break;
                }
            }
            if (include == false)
            {
                break;
            }
        }

        if (include)
        {
            proximal_targets[query] = proximals;
        }
    }



    vector<Crispr> all_crisprs;

    for ( auto const& [query, proximal] : proximal_targets)
    {
        vector<Crispr> candidates;
        for (ui k = Prospector::k_end-1; k >= Prospector::k_start; k--)
        {
            ui allowed_mutations = k / Prospector::mutant_tolerance_ratio;

            vector<ui> genome_indices;

            genome_indices.push_back(query);

            for (ui target : proximal)
            {
                ui end = genome_indices[genome_indices.size()-1] + k;
                if (target < end || target - end < Prospector::spacer_min) continue; // || guards against overflow
                if (target - end > Prospector::spacer_max) break;

                if (mutant(genome.c_str(), encoding.encoding, k, allowed_mutations, query, target))
                    genome_indices.push_back(target);

            }

            if (genome_indices.size() >= Prospector::repeats_min)
            {
                Crispr c(k, genome_indices, genome_indices.size());
                candidates.push_back(c);
            }
        }
   
        if (candidates.size() > 0)
        {
            Util::sort(candidates, [](const Crispr& a, const Crispr& b) {  return a.size > b.size; } );
            ui best_size = candidates[0].size;
            for (const Crispr& crispr : candidates)
            {
                if (crispr.size != best_size)
                    break;

                all_crisprs.push_back(crispr);
            }
        }
    }

    time(start, "post-kernel crispr generation");
    fmt::print("\tprospector returned {} crisprs\n", all_crisprs.size());
    return all_crisprs;
}


vector<Crispr> get_crisprs(const string& genome)
{
    vector<Crispr> crisprs = prospector_main(genome);      
    CrisprUtil::cache_crispr_information(genome, crisprs);
    // CrisprUtil::debug(crisprs, genome, 1085431-1000, 1086854+1000);
    crisprs = Util::filter(crisprs, [](const Crispr& c) { return c.overall_heuristic >= 0.75; });
    Util::sort(crisprs, CrisprUtil::heuristic_greater);
    crisprs = CrisprUtil::get_domain_best(crisprs);
    Util::sort(crisprs, [](const Crispr& a, const Crispr&b) { return a.start < b.start; });
    return crisprs;
}


void stdrun(const string& genome, const vector<CasProfile>& cas_profiles)
{
    auto start = time();

    vector<Crispr> crisprs = get_crisprs(genome);
    vector<Translation> flanks = CasUtil::get_translations(genome, crisprs);
    vector<Fragment> fragments = CasUtil::cas(cas_profiles, flanks, genome);

    CrisprUtil::print(genome, crisprs);
    CasUtil::print_fragments(crisprs, fragments, genome);

    time(start, "stdrun");
    fmt::print("\n\n");
}

void write(string cas_file, string cache_file)
{
    vector<CasProfile> cas_profiles = CasUtil::load(cas_file, CasUtil::gen_n);
    CasUtil::write_cache(cache_file, cas_profiles);
    fmt::print("cache written, exiting...\n");
    exit(0);
}

vector<CasProfile> read(string cas_file, string cache_file)
{
    CasUtil::load_cache("crispr-data/cas/cache.fasta");
    vector<CasProfile> cas_profiles = CasUtil::load(cas_file, CasUtil::get_n);
    return cas_profiles;
}

int main()
{
    Prospector::device_init();

    printf("running main...\n");
    auto start = time();

    string genome_dir = "crispr-data/genome";
    string cas_file = "crispr-data/cas/cas.fasta";
    string cache_file = "crispr-data/cas/cache.fasta";

    map<string, string> genomes = Util::load_genomes(genome_dir);

    // write(cas_file, cache_file);
    vector<CasProfile> cas_profiles = read(cas_file, cache_file);

    // stdrun(genomes["thermophilus"], cas_profiles);
    stdrun(genomes["pyogenes"], cas_profiles);
    time(start, "main");

    return 0;                                                                                                           
}

