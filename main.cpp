//std
#include "stdafx.h"
#include "fmt/core.h"
#include "fmt/format.h"
#include <bitset>
#include <cstddef>
#include <filesystem>
namespace fs = std::filesystem;

//proj
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"
#include "cas.h"

void finish()
{
    fmt::print("terminate called\n");
    exit(0);
}

void debug(vector<Crispr> crisprs, string genome, ui start, ui end)
{

    vector<Crispr> filtered = filter(crisprs, [&](Crispr c) { return c.start > start-100 && c.end < end+100; } );

    sort(filtered.begin(), filtered.end(), CrisprUtil::heuristic_less);

    int how_many = filtered.size();
    for (size_t i = filtered.size()-how_many; i < filtered.size(); i++)
        filtered[i].print(genome);
    finish();
}

vector<string> load_genomes(string dir)
{
    vector<string> genomes;
    for (const auto& entry : filesystem::directory_iterator(dir))
        genomes.push_back(parse_fasta_single(entry.path()));
    return genomes;
}


map<char, ull> _scheme {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}
};
const ull bits = 2;

ull compute(ull e1, ull e2)
{
    ull _xor = (e1 ^ e2);
    ull evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ull oddBits = _xor & 0x5555555555555555ull;
    ull comp = (evenBits >> 1) | oddBits;
    return popcount(comp);
}


ull encode(string kmer)
{
    ull encoding = 0;
    for (int i = 0; i < kmer.size(); i++)
        encoding += _scheme[kmer[i]] << i * bits;
    return encoding;
}

void print_encoding(ull encoding)
{
    bitset<64> a(encoding);
    cout << a << endl;
}

void stdrun(const string& genome, string cas_dir)
{
    double start = omp_get_wtime();

    vector<Crispr> crisprs = Prospector::prospector_main(genome);      

    CrisprUtil::cache_crispr_information(genome, crisprs);
    vector<Crispr> good = filter(crisprs, [](const Crispr& c) { return c.overall_heuristic >= 0.75; });
    sort(good.begin(), good.end(), CrisprUtil::heuristic_greater);
    vector<Crispr> final = CrisprUtil::get_domain_best(good);
    sort(final.begin(), final.end(), [](const Crispr& a, const Crispr&b) { return a.start < b.start; });
    CrisprUtil::print(genome, final);


    double _start = omp_get_wtime();

    vector<Translation> downstreams;
    vector<Translation> upstreams;

    for (const Crispr& c : final)
    {
        downstreams.push_back(Translation::from_crispr_down(genome, c));
        upstreams.push_back(Translation::from_crispr_up(genome, c));
    }

    done(_start, "translation gen");

    vector<Fragment> fragments = Cas::cas(genome, final, cas_dir, downstreams, upstreams);
    Cas::print_fragments(final, fragments);
    done(start, "stdrun");
}

int main()
{
    

    printf("running main...\n");
    double start = omp_get_wtime();

    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";
    vector<string> genomes = load_genomes(genome_dir);

    stdrun(genomes[0], cas_dir);
    stdrun(genomes[1], cas_dir);


    done(start, "main");

    finish();
    return 0;                                                                                                           
}

