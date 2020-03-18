#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"


int main()
{
    printf("running invoker...\n");
    double start = omp_get_wtime();

    map<string, string> genomes =
    {
            {"thermophilus", parse_fasta_single("crispr-data/streptococcus_thermophilus.fasta")},
            {"pyogenes", parse_fasta_single("crispr-data/pyogenes.fasta")}
    };

    string genome = genomes["thermophilus"];

    vector<Crispr> crisprs = Prospector::prospector_main(genome);
    
    CrisprUtil::cache_crispr_information(crisprs, genome);
    
    vector<Crispr> good_heuristic_crisprs;
    for (const Crispr& crispr : crisprs)
    {
        if (crispr.overall_heuristic > 0.5)
            good_heuristic_crisprs.push_back(crispr);
    }

    sort(good_heuristic_crisprs.begin(), good_heuristic_crisprs.end(), CrisprUtil::heuristic_greater);
    
    vector<Crispr> domain_best = CrisprUtil::get_domain_best(good_heuristic_crisprs);

    map<string, int> spacer_scores = CrisprUtil::get_spacer_scores(domain_best);
    
    vector<Crispr> final = CrisprUtil::spacer_score_filtered(domain_best, spacer_scores);

    CrisprUtil::print(genome, final, spacer_scores);

    CrisprUtil::cas(genome, final, 5, 10000);
    
    done(start, "invoker");

    return 0;
}

