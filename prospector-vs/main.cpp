#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

//string crispr_type(vector<Gene*>& genes)
//{
//    static const map<string, string> class_lookup
//    {
//        {"I",   "1"},
//        {"III", "1"},
//        {"IV",  "1"},
//        {"II",  "2"},
//        {"V",   "2"},
//        {"VI",  "2"},
//        {"?",   "?"},
//    };
//
//    static const map<string, string> type_lookup
//    {
//        {"cas3",  "I"},
//        {"cas10", "III"},
//        {"cas8",  "IV"},
//        {"cas9",  "II"},
//        {"cas12", "V"},
//        {"cas13", "VI"},
//    };
//
//    for (Gene* gene : genes)
//    {
//        if (type_lookup.contains(gene->gn))
//        {
//            return type_lookup.at(gene->gn);
//        }
//    }
//    return "?";
//}

static vector<CasProfile*> profiles;

string Fragment::to_string_debug()
{
    string amino_family = Util::translate_genome(*reference_genome, genome_begin, genome_final, reference_translation->pos);
    string amino_cds = Util::translate_genome(*reference_genome, expanded_genome_begin, expanded_genome_final, reference_translation->pos);

    string dna_family = reference_genome->substr(genome_begin, genome_final - genome_begin);
    string dna_cds = reference_genome->substr(expanded_genome_begin, expanded_genome_final - expanded_genome_begin);

    std::ostringstream out;
    out << fmt::format("{}\n", reference_profile->gn);
    out << fmt::format("\t{}\n", reference_translation->pos ? "+" : "-");
    out << fmt::format("\t{}\n", reference_translation->reference_crispr->identifier_string());
    out << fmt::format("\t{}...{}\n", genome_begin, genome_final);
    out << fmt::format("\t{}...{}\n", expanded_genome_begin, expanded_genome_final);
    out << fmt::format("\t{}\n", amino_family);
    out << fmt::format("\t{}\n", amino_cds);
    out << fmt::format("\t{}\n", dna_family);
    out << fmt::format("\t{}\n", dna_cds);
    return out.str();
}

string Fragment::to_string_summary()
{
    return fmt::format("{}\t{}\t{}\t{}\n", expanded_genome_begin + 1, expanded_genome_final, reference_translation->pos ? "+" : "-", reference_profile->gn);
}

string Crispr::to_string_debug()
{
    std::ostringstream out;

    out << fmt::format("{} - {} {}\n", start, end, k);
    out << fmt::format("{}h {}r {}s {}v\n", overall_heuristic, conservation_repeats, conservation_spacers2, spacer_variance);
    out << fmt::format("\t{} repeats\n", repeats.size());
    out << fmt::format("\t{} spacers\n", spacers.size());

    for (ull i = 0; i < repeats.size(); i++)
    {
        string repeat = repeats[i];
        int mismatches = Util::mismatch_count(repeat);
        int matches = repeat.length() / 2 - mismatches;
        double score = (double)matches / (double)(repeat.length() / 2);
        int start = genome_indices[i];
        int end = start + k - 1;
        int dist = i == 0 ? 0 : genome_indices[i] - (genome_indices[i - 1] + k);

        //printf("\t\t");
        //printf("%d/%zd", matches, repeat.length()/2);
        //printf(" %d %s %d", start, repeat.c_str(), end);
        //printf(" %d", dist);
        //printf(" %f", score);

        out << fmt::format("\t\t{} {} {} {}\n", start, repeat, end, i < repeats.size() - 1 ? spacers[i] : "");
    }

    out << "\n";

    return out.str();
}

string Crispr::to_string_summary()
{
    return fmt::format("{}\t{}\t{}\t{}\t{}h\n", start + 1, end + 1, "?", "CRISPR", overall_heuristic);
}

void write_all(const string& genome_name, const string& genome, const vector<Crispr*>& crisprs, vector<Fragment*>& fragments, std::ofstream& results)
{
    results << fmt::format("===\t{}\n", genome_name);

    vector<Locus*> loci;
    for (Crispr* c : crisprs) loci.push_back(c);
    for (Fragment* f : fragments) loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) {return a->get_start() < b->get_start(); });

    for (Locus* l : loci)
        results << l->to_string_summary();
}

void prospect_genome(string genome_path, std::ofstream& results)
{
    fmt::print("\n\n\n");

    auto start = time();

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    write_all(genome_path, genome, crisprs, fragments, results);

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
    // ------ start -----------
    auto start_main = time();


    // ---- serialziation ------
    //CasProfileUtil::pfam_filter("T:\\data\\Pfam-A.seed", "T:\\data\\Pfam-A.filt");
    //vector<const CasProfile*> profiles = CasProfileUtil::pfam("T:\\crispr\\cas\\Pfam-A.full_filt");
    //CasProfileUtil::serialize("T:\\crispr\\cas\\serial_staging", profiles);
    //exit(0);

    // -------- init -----------
    Prospector::device_init();
    profiles = CasProfileUtil::deserialize("T:\\crispr\\cas\\serial");

    // ----------- debug --------
    //serialization();
    //Debug::translation_print("T:\\crispr\\supp\\genomes\\GCA_000011125.1_ASM1112v1_genomic.fna", 784236, 784515, false, 10);
    //exit(0);

    // ----------- results --------------
    std::ofstream results("results.txt");
    prospect_genome_dir("T:\\crispr\\supp\\genomes", results);
    //prospect_genome("T:\\crispr\\supp\\genomes\\GCA_000011125.1_ASM1112v1_genomic.fna", results);
    results.close();

    // --------- close -----------
    start_main = time(start_main, "main"); return 0;                                                                                                           
}


