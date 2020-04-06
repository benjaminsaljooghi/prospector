//std
#include "stdafx.h"

//proj
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"
#include "cas.h"


map<char, ui> scheme {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3} 
};


ui encoded(const string& kmer)
{
    #if DEBUG == 1
    assert(kmer.size() == SIZE);
    #endif
    ui e = 0;
    for (int i = 0; i < kmer.size(); i++)
        e |= scheme.at(kmer[i]) << (i * BITS);
    return e;
}

ui* encoded_genome(const string& genome)
{
    double __start = omp_get_wtime();
    ui num = genome.size() - SIZE + 1;
    ui* encoding = (ui*) malloc(sizeof(ui) * num);
    #pragma omp parallel for
    for (ui i = 0; i < num; i++) encoding[i] = encoded(genome.substr(i, SIZE));
    done(__start, "genome encoding");
    return encoding;
}

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
    double start = omp_get_wtime();
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
    done(start, "q_substrate");
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



vector<vector<ui>> single_k_from_q_substrate(const char* genome, const vector<ui>& queries, ui* genome_encoding, const ui& k)
{
    vector<vector<ui>> crisprs;
    ui allowed_mutations = k / MUTANT_TOLERANCE_RATIO;

    for (ui _q = 0; _q < queries.size(); _q++)
    {
        ui q = queries[_q];

        vector<ui> crispr;
        crispr.push_back(q);

        ui bound = q + k + SPACER_SKIP;
        
        for (ui t = bound; t - bound <= SPACER_MAX; t++)
        {
            if (mutant(genome, genome_encoding, k, allowed_mutations, q, t))
            {
                crispr.push_back(t);
                bound = t + k + SPACER_SKIP;
                t = bound;
            }
        }
        crisprs.push_back(crispr);
    }
    return crisprs;
}



vector<Crispr> prospector_main(const string& genome)
{
    ui* genome_encoding = encoded_genome(genome);
    ui genome_encoding_size = genome.size() - SIZE + 1;

    unsigned char* qmap = Prospector::get_qmap(genome_encoding, genome_encoding_size);

    vector<ui> queries = q_substrate(qmap, genome_encoding_size);
    
    double total_single_k_time = 0;
    vector<Crispr> all_crisprs;
    for (ui k = K_START; k < K_END; k++)
    {
        double __start = omp_get_wtime();
        vector<vector<ui>> crisprs = single_k_from_q_substrate(genome.c_str(), queries, genome_encoding, k);
        double __end = omp_get_wtime();
        total_single_k_time += __end - __start;

        for (vector<ui> c : crisprs)
        {
            if (c.size() >= MIN_REPEATS)
            {
                Crispr _c(k, c, c.size());
                all_crisprs.push_back(_c);   
            }
        }
    }
    fmt::print("\tprospector returned {} crisprs\n", all_crisprs.size());
    return all_crisprs;
}




vector<Crispr> get_crisprs(const string& genome)
{
    vector<Crispr> crisprs = prospector_main(genome);      

    CrisprUtil::cache_crispr_information(genome, crisprs);

    crisprs = filter(crisprs, [](const Crispr& c) { return c.overall_heuristic >= 0.75; });

    sort(crisprs, CrisprUtil::heuristic_greater);

    crisprs = CrisprUtil::get_domain_best(crisprs);

    sort(crisprs, [](const Crispr& a, const Crispr&b) { return a.start < b.start; });

    return crisprs;
}

void stdrun(const string& genome, const vector<CasProfile>& cas_profiles)
{
    double start = omp_get_wtime();

    vector<Crispr> crisprs = get_crisprs(genome);

    vector<Translation> downstreams;
    vector<Translation> upstreams;

    for (const Crispr& c : crisprs) downstreams.push_back(Translation::from_crispr_down(genome, c));
    for (const Crispr& c: crisprs) upstreams.push_back(Translation::from_crispr_up(genome, c));
    
    vector<Fragment> fragments = Cas::cas(genome, crisprs, cas_profiles, downstreams, upstreams);

    // CrisprUtil::print(genome, crisprs);
    // Cas::print_fragments(crisprs, fragments);

    done(start, "stdrun");
    fmt::print("\n\n");
}

int main()
{
    printf("running main...\n");
    double start = omp_get_wtime();

    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";
    vector<string> genomes = Util::load_genomes(genome_dir);
    vector<CasProfile> cas_profiles = CasProfile::load_casprofiles(cas_dir, K_FRAGMENT);

    stdrun(genomes[0], cas_profiles);
    stdrun(genomes[1], cas_profiles);

    done(start, "main");

    return 0;                                                                                                           
}

