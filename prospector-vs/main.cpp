#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

void prospect(const vector<const CasProfile*>& cas_profiles, map<string, string>& genomes, const string& genome_name)
{
    auto start_run = time();
    string& genome = genomes.at(genome_name);
    
    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(cas_profiles, translations, genome);
    map<string, vector<Gene*>> genes = Cas::assemble_genes(crisprs, fragments);
    Cas::print_all(crisprs, genes, genome);
    start_run = time(start_run, genome_name.c_str());
}

string genome_dir = "T:\\crispr\\genome";
string serial_dir = "T:\\crispr\\cas\\serial";

void instantiation_routine()
{
    //CasProfileUtil::pfam_filter("T:\\data\\Pfam-A.seed", "T:\\data\\Pfam-A.filt");
    vector<const CasProfile*> profiles = CasProfileUtil::pfam("T:\\crispr\\cas\\Pfam-A.filt");
    CasProfileUtil::serialize(serial_dir, profiles);
    exit(0);
}

int main()
{
    Prospector::device_init(); 
    auto profiles = CasProfileUtil::deserialize(serial_dir);
    auto genomes = Util::load_genomes(genome_dir);

    auto start_main = time();

    prospect(profiles, genomes, "thermosaccharolyticum");

    start_main = time(start_main, "main");
    return 0;                                                                                                           
}


