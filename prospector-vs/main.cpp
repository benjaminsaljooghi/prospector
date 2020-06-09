#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

void serialization()
{
    CasProfileUtil::pfam_filter("T:\\data\\Pfam-A.seed", "T:\\data\\Pfam-A.filt");
    vector<const CasProfile*> profiles = CasProfileUtil::pfam("T:\\crispr\\cas\\Pfam-A.filt");
    CasProfileUtil::serialize("T:\\crispr\\cas\\serial", profiles);
    exit(0);
}

vector<const CasProfile*> profiles;

void init()
{
    Prospector::device_init();
    profiles = CasProfileUtil::deserialize("T:\\crispr\\cas\\serial");
}

void prospect_genome(string genome_path)
{
    string genome = Util::load_genome(genome_path);
    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    map<string, vector<Gene*>> genes = Cas::assemble_genes(crisprs, fragments);
    Cas::print_all(crisprs, genes, genome);
}

void prospect_genome_dir(string genome_dir)
{
    ui i = 0;
    for (const auto& entry : filesystem::directory_iterator(genome_dir))
    {
        if (i++ > 10) break;
        string genome_path = entry.path().string();
        prospect_genome(genome_path);
    }
}

int main()
{
    init();

    auto start_main = time();
    prospect_genome_dir("T:\\crispr\\supp\\genomes");
    start_main = time(start_main, "main");

    return 0;                                                                                                           
}


