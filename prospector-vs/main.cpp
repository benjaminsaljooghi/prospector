#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"


vector<CasProfile*> profiles;


map<string, string> class_lookup
{
    {"I",   "1"},
    {"III", "1"},
    {"IV",  "1"},
    {"II",  "2"},
    {"V",   "2"},
    {"VI",  "2"},
    {"?",   "?"},
};

map<string, string> type_lookup
{
    {"cas3",  "I"},
    {"cas10", "III"},
    {"cas8",  "IV"},
    {"cas9",  "II"},
    {"cas12", "V"},
    {"cas13", "VI"},
};

//string crispr_type(vector<Gene*>& genes)
//{
//    for (Gene* gene : genes)
//    {
//        if (type_lookup.contains(gene->gn))
//        {
//            return type_lookup.at(gene->gn);
//        }
//    }
//    return "?";
//}

//void print_gene_summary(Gene& gene)
//{
//    fmt::print("\t{}\n", gene.reference_profile->gn);
//    for (const Fragment& f : gene.fragments)
//    {
//        fmt::print("\t\t{}...{}\n", f.details->genome_start, f.details->genome_final);
//    }
//    fmt::print("\n");
//}

//void print_fragment_debug(const Fragment* f, const string& genome)
//{
//

//}

//void print_gene_debug(Gene* gene, const string& genome)
//{
//    //ui size = gene->size();
//    fmt::print("\t{}:{}\n", gene->gn, gene->fragments[0]->reference_translation->pos);
//    for (const Fragment* fragment : gene->fragments)
//    {
//        print_fragment_debug(fragment, genome);
//    }
//
//    fmt::print("\n");
//}


void serialization()
{
    CasProfileUtil::pfam_filter("T:\\data\\Pfam-A.seed", "T:\\data\\Pfam-A.filt");
    vector<const CasProfile*> profiles = CasProfileUtil::pfam("T:\\crispr\\cas\\Pfam-A.filt");
    CasProfileUtil::serialize("T:\\crispr\\cas\\serial", profiles);
    exit(0);
}


void init()
{
    Prospector::device_init();
    profiles = CasProfileUtil::deserialize("T:\\crispr\\cas\\serial");
}

void print_all(const string& genome_name, const string& genome, const vector<Crispr*>& crisprs, vector<Fragment*>& fragments, std::ofstream& results)
{

    results << fmt::format("===\t{}\n", genome_name);

    std::sort(fragments.begin(), fragments.end(), [](Fragment* a, Fragment* b) {return a->expanded_genome_begin < b->expanded_genome_begin; });


    auto debug = [&genome](Fragment* f) {
        string amino_family = Util::translate_genome(genome, f->genome_begin, f->genome_final, f->reference_translation->pos);
        string amino_cds = Util::translate_genome(genome, f->expanded_genome_begin, f->expanded_genome_final, f->reference_translation->pos);

        string dna_family = genome.substr(f->genome_begin, f->genome_final - f->genome_begin);
        string dna_cds = genome.substr(f->expanded_genome_begin, f->expanded_genome_final - f->expanded_genome_begin);

        std::ostringstream out;

        out << fmt::format("{}\n", f->reference_profile->gn);

        out << fmt::format("\t{}\n", f->reference_translation->reference_crispr->identifier_string());

        out << fmt::format("\t{}...{}\n", f->genome_begin, f->genome_final);
        out << fmt::format("\t{}...{}\n", f->expanded_genome_begin, f->expanded_genome_final);

        out << fmt::format("\t{}\n", amino_family);
        out << fmt::format("\t{}\n", amino_cds);

        out << fmt::format("\t{}\n", dna_family);
        out << fmt::format("\t{}\n", dna_cds);

        return out.str();
    };

    auto summary = [](Fragment* f) {
        return fmt::format("{}\t{}\t{}\t{}\n", f->expanded_genome_begin, f->expanded_genome_final, f->reference_translation->pos ? "+" : "-", f->reference_profile->gn);
    };

    auto c_summary = [](Crispr* c) {
        return fmt::format("{}\t{}\n", c->start, c->end);
    };

    for (Crispr* c : crisprs)
        results << c_summary(c);

    for (Fragment* f : fragments)
        results << summary(f);

}

void prospect_genome(string genome_path, std::ofstream& results)
{
    fmt::print("\n\n\n");

    auto start = time();

    string genome = Util::load_genome(genome_path);

    //Debug::translation_print(genome, 776022, 776754, false, 10);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    print_all(genome_path, genome, crisprs, fragments, results);


    start = time(start, genome_path.c_str());
}

void prospect_genome_dir(string genome_dir, std::ofstream& results)
{

    ui i = 0;
    for (const auto& entry : filesystem::directory_iterator(genome_dir))
    {
        if (i++ > 10) break;
        string genome_path = entry.path().string();
        prospect_genome(genome_path, results);
    }

}

int main()
{
    init();

    auto start_main = time();

    std::ofstream results("results.txt");

    prospect_genome_dir("T:\\crispr\\supp\\genomes", results);
    //prospect_genome("T:\\crispr\\supp\\genomes\\GCA_000011125.1_ASM1112v1_genomic.fna", results);

    results.close();

    start_main = time(start_main, "main");

    return 0;                                                                                                           
}


