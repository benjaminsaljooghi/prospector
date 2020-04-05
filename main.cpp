//std
#include "stdafx.h"



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


void stdrun(const string& genome, string cas_dir)
{
    double start = omp_get_wtime();

    vector<Crispr> crisprs;
    // crisprs = Prospector::prospector_main(genome);      

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
    vector<string> genomes = Util::load_genomes(genome_dir);

    stdrun(genomes[0], cas_dir);
    stdrun(genomes[1], cas_dir);


    done(start, "main");

    finish();
    return 0;                                                                                                           
}

