#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"

map<string, int> get_spacer_scores(vector<Crispr>& crisprs)
{
    set<string> all_spacers;
    for (Crispr& crispr : crisprs)
        all_spacers.insert(crispr.spacers.begin(), crispr.spacers.end());
    return BLAST(all_spacers);
}


int main()
{
    printf("running invoker...\n");
    double start = omp_get_wtime();

    map<string, string> genomes = {
            {"thermophilus", parse_fasta("crispr-data/streptococcus_thermophilus.fasta").begin()->second},
            {"pyogenes", parse_fasta("crispr-data/pyogenes.fasta").begin()->second}
    };
    string genome = genomes["thermophilus"];

    vector<Crispr> crisprs = Prospector::prospector_main(genome);
    CrisprUtil::cache_crispr_information(crisprs, genome);
    sort(crisprs.begin(), crisprs.end());
    vector<Crispr> domain_best = CrisprUtil::get_domain_best(crisprs);
    map<string, int> spacer_scores = get_spacer_scores(domain_best);
    vector<Crispr> final = CrisprUtil::score_filtered(domain_best, spacer_scores);
    CrisprUtil::print(genome, final, spacer_scores);
    cas(genome, final, 5, 10000);
    done(start, "invoker");
    return 0;
}

