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


vector<ui> q_substrate(unsigned char* qmap, ui genome_encoding_size)
{
    // how many qs in this map are containment oriented
    auto start = time();
    // ui count = 0;
    vector<ui> queries;
    for (ui query = 0; query < genome_encoding_size - 200; query++)
    {
        for (ui i = 0; i < MAP_SIZE; i++)
        {
            if (qmap[(query*MAP_SIZE) + i] <= (SIZE / MUTANT_TOLERANCE_RATIO)) // 1 because 8 / 5 = 1
            {
                queries.push_back(query);
                break;
            }
        }
    }
    time(start, "q_substrate");
    // printf("%d %zd\n", genome_encoding_size-200, queries.size());
    // return count;
    return queries;
}



bool mutant(const char* genome, const ui* genome_encoding, const ui& k, const ui& allowed_mutations, const ui& i, const ui& j)
{
    ui diff = 0;
    const ui chunks = k / SIZE;
    // may generate a lot of crisprs that are filtered later (expensive) given that SIZE is large (16) here.
    // option is to do a more accurate mutation calculation either using per-char post the chunk division
    // or to encode entire kmers up to MAP_SIZ into ull's to compute the difference efficiently.
    // post k=MAP_SIZ we can use the MAP_SIZ-long ull initially, and then compute a per-char difference afterwards.

    for (ui chunk = 0; chunk < chunks; chunk++)
    {
        ui _i = genome_encoding[i + (chunk * SIZE)];
        ui _j = genome_encoding[j + (chunk * SIZE)];
        diff += difference_cpu(_i, _j);
        if (diff > allowed_mutations)
        {
            return false;
        }
    }
    const ui checked_so_far = (chunks * SIZE);

    return diff <= checked_so_far / MUTANT_TOLERANCE_RATIO;

    // for (ui __i = checked_so_far; i < k; __i++)
    // {
        // diff += genome[i + checked_so_far + __i] == genome[j + checked_so_far + __i] ? 0 : 1; 
    // }
    // return diff <= allowed_mutations;
    
}


void debug_map()
{
    // ui query = 1283501;
    // ui q = genome_encoding[query];
    // for (ui i = 0; i < 1000; i++)
    // {
    //     ui pos = query + K_START + SPACER_SKIP + i;
    //     ui diff = difference_cpu(genome_encoding[query], genome_encoding[pos]);

    //     printf("%s %d %d\n", genome.substr(pos, SIZE).c_str(), pos, diff);
    // }
}



// vector<vector<ui>> single_k_from_q_substrate(const char* genome, const vector<ui>& queries, ui* genome_encoding, const ui& k)
// {
//     vector<vector<ui>> crisprs;
//     ui allowed_mutations = k / MUTANT_TOLERANCE_RATIO;

   
//     for (ui q : queries)
//     {
//         // is this query contained within the previous crispr?
//         // bool cont = false;
//         // for (vector<ui> c : crisprs)
//         // {
//         //     if (contains(c, q))
//         //     {
//         //         cont = true;
//         //         break;
//         //     }
//         // }
//         // if (cont) continue;

//         vector<ui> crispr;
//         crispr.push_back(q);

//         ui bound = q + k + SPACER_SKIP;
        
//         for (ui t = bound; t - bound <= SPACER_MAX; t++)
//         {
//             if (mutant(genome, genome_encoding, k, allowed_mutations, q, t))
//             {
//                 crispr.push_back(t);
//                 bound = t + k + SPACER_SKIP;
//                 t = bound;
//             }
//         }
//         crisprs.push_back(crispr);
//     }
//     return crisprs;
// }



vector<Crispr> prospector_main(const string& genome)
{
    Prospector::Encoding encoding = Prospector::get_genome_encoding(genome.c_str(), genome.size());

    uc* qmap = Prospector::get_qmap(encoding.d_encoding, encoding.size);

    vector<ui> queries = q_substrate(qmap, encoding.size);
    
    uc* qmap_3000 = Prospector::get_qmap3000(encoding.d_encoding, encoding.size, &queries[0], queries.size());


    // per q_i proximal targets

    vector<vector<ui>> proximal_targets;
    ui tolerance = 16 / MUTANT_TOLERANCE_RATIO;

    for (ui q_i = 0; q_i < queries.size(); q_i++)
    {
        vector<ui> proximals;
        for (ui i = 0; i < 3000; i++)
        {
            if (qmap_3000[q_i * 3000 + i] <= tolerance)
            {
                proximals.push_back(i);
            }
        }
        proximal_targets.push_back(proximals);
    }


    auto start = time();


    vector<Crispr> all_crisprs;

    for (ui k = K_START; k < K_END; k++)
    {
        // build the largest crispr, and then leap to the next crispr
        for (ui query_i = 0; query_i < queries.size(); query_i++)
        {
            ui query = queries[query_i];

            // form largest crispr by iterating over the targets of interest
            

        }

    }

    // vector<Crispr> all_crisprs;
    // for (ui k = K_START; k < K_END; k++)
    // {
    //     vector<vector<ui>> crisprs = single_k_from_q_substrate(genome.c_str(), queries, encoding.encoding, k);

    //     for (vector<ui> c : crisprs)
    //     {
    //         if (c.size() >= MIN_REPEATS)
    //         {
    //             Crispr _c(k, c, c.size());
    //             all_crisprs.push_back(_c);   
    //         }
    //     }
    // }
    time(start, "crisps from q_substrate");

    fmt::print("\tprospector returned {} crisprs\n", all_crisprs.size());
    return all_crisprs;
}





vector<Crispr> get_crisprs(const string& genome)
{
    vector<Crispr> crisprs = prospector_main(genome);      

    // how many crisprs are a perfect subset in that they contain an equal k value
    // but just fewer repeats than another crispr in that all other crispr genome indices are the same
    
    // auto start = time();
    // ui count = 0;
    // vector<Crispr> __subset;
    // for (ui i = 0; i < crisprs.size(); i++)
    // {
    //     bool include = true;
    //     const Crispr& a = crisprs[i];
    //     for (ui j = 0; j < crisprs.size(); j++)
    //     {
    //         const Crispr& b = crisprs[j];

    //         if (i == j) continue;
    //         if (a.k != b.k) continue;

    //         // check that all starts in a are contained within b
    //         if (subset(a.genome_indices, b.genome_indices))
    //         {
    //             // subset.push_back(a);
    //             include = false;
    //             break;
    //         }
    //     }
    //     if (include)
    //     {
    //         __subset.push_back(a);
    //     }
    // }
    // time(start, "subset removal");

    // fmt::print("{} {}\n", __subset.size(), crisprs.size());
    // crisprs = __subset;


    CrisprUtil::cache_crispr_information(genome, crisprs);


    // CrisprUtil::debug(crisprs, genome, 1824943, 1827551);


    crisprs = filter(crisprs, [](const Crispr& c) { return c.overall_heuristic >= 0.75; });

    sort(crisprs, CrisprUtil::heuristic_greater);

    crisprs = CrisprUtil::get_domain_best(crisprs);

    sort(crisprs, [](const Crispr& a, const Crispr&b) { return a.start < b.start; });

    return crisprs;
}


void stdrun(const string& genome, const vector<CasProfile>& cas_profiles)
{
    auto start = time();

    vector<Crispr> crisprs = get_crisprs(genome);
    vector<Flanks> flanks = CasUtil::get_flanks(genome, crisprs);
    vector<Fragment> fragments = CasUtil::cas(genome, crisprs, cas_profiles, flanks);

    CrisprUtil::print(genome, crisprs);
    CasUtil::print_fragments(crisprs, fragments);

    time(start, "stdrun");
    fmt::print("\n\n");
}

int main()
{
    Prospector::device_init();


    printf("running main...\n");
    auto start = time();

    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";
    vector<string> genomes = Util::load_genomes(genome_dir);
    vector<CasProfile> cas_profiles = CasUtil::load(cas_dir, K_FRAGMENT);

    // stdrun(genomes[0], cas_profiles);
    stdrun(genomes[1], cas_profiles);

    time(start, "main");

    return 0;                                                                                                           
}

