#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"
#include "cas_profiles.h"
#include "array_discovery.h"

#define DEBUG 1

static map<string, string> domain_to_gn;
string genome_dir = "T:\\crispr\\genome";
string gn_path = "T:\\crispr\\groundtruth\\gn.tsv";
string serial = "T:\\crispr\\cas\\serial";
static vector<CasProfile*> profiles;
static map<string, CasProfile*> profiles_organized;
std::ofstream out_debug("out_debug.txt");
std::ofstream out_gene("out_gene.txt");
std::ofstream out_domain("out_domain.txt");

void init()
{
    Prospector::device_init();
    profiles = CasProfileUtil::deserialize(serial);

    for (CasProfile* profile : profiles)
    {
        profiles_organized[profile->gn] = profile;
    }

    std::ifstream a(gn_path);
    string line;
    while (getline(a, line))
    {
        auto split = Util::parse(line, "\t");
        domain_to_gn[split[0]] = split[1];
    }


}

bool Fragment::is_domain()
{
    return !domain_to_gn.contains(this->reference_profile->gn);
}

bool Fragment::is_gene()
{
    return !this->is_domain();
}

void sage_interpreter(string path, string genome_dir)
{
    std::ifstream in(path);
    string line;
    string genome_accession = "";
    string genome = "";
    std::ofstream sage_interpretation("sage_interpretation.txt");

    while (std::getline(in, line))
    {
        if (line[0] == '-')
        {
            continue;
        }

        auto split = Util::parse(line, "\t");

        if (split[0] == "===")
        {
            genome_accession = split[1];

            for (const auto& entry : filesystem::directory_iterator(genome_dir))
            {
                string genome_path = entry.path().string();
            
                if (genome_path.find(genome_accession) != string::npos)
                {
                    genome = Util::load_genome(genome_path);
                    break;
                }
            }

            continue;
        }
        

        auto gen_debug_str = [&genome](string kind, ui begin, ui final, string strand) {
            auto a = genome.substr(begin, final - begin);

            auto b = kind == "CRISPR" ? "" : Debug::translation_test(genome, begin, final, strand == "+", 0);
        
            std::ostringstream stream;

            stream << fmt::format("\t{}\t{}..{}\t{}\n", kind, begin, final, strand);
            stream << "\t" << a << "\n";
            stream << "\t" << b << "\n";

            return stream.str();
        };
     

        if (split[0] != "")
        {
            ui g_begin = stoi(split[0]);
            ui g_final = stoi(split[1]);
            string g_strand = split[2];
            string kind = split[3];
            string str = gen_debug_str(kind, g_begin, g_final, g_strand);
            sage_interpretation << fmt::format("ground:\n{}\n", str);
        }
        else
        {
            sage_interpretation << "ground:\nno data\n";
        }

        if (split[5] != "")
        {
            ui t_begin = stoi(split[5]);
            ui t_final = stoi(split[6]);
            string t_strand = split[7];
            string kind = split[8];
            auto str = gen_debug_str(kind, t_begin, t_final, t_strand);
            sage_interpretation << fmt::format("target:\n{}\n", str);
        }
        else
        {
            sage_interpretation << "target:\nno data\n";
        }

        sage_interpretation << "\n\n";
    }

    sage_interpretation.close();
}

string get_accession(string genome_path)
{
    static const std::regex accession_pattern("GC[AF]_[0-9]+\.[0-9]+", regex_constants::icase);
    std::filesystem::path my_path(genome_path);
    string stem = my_path.stem().string();
    std::cmatch match;
    bool result = std::regex_search(stem.c_str(), match, accession_pattern);
    return match.str();
}

string get_header(string genome_path)
{
    return fmt::format("===\t{}\t\t{}\n", genome_path, get_accession(genome_path));
}

void write_loci(string genome_path, vector<Locus*> loci)
{
    string header = get_header(genome_path);

    out_gene << header;
    out_debug << header;
    out_domain << header;

    for (Locus* l : loci)
    {
        string debug_info = l->to_string_debug();
        string summary_info = l->to_string_summary();

        out_debug << debug_info;

        if (l->is_crispr())
        {
            out_domain << summary_info;
            out_gene << summary_info;
        }

        if (l->is_domain())
        {
            out_domain << summary_info;
        }
        
        if (l->is_gene())
        {
            out_gene << summary_info;
        }
    }
}

void prospect_genome(string genome_path)
{
    fmt::print("\n\n\n");

    auto start = time();
    
    string genome = Util::load_genome(genome_path);
    start = time(start, "load genome");

    vector<Crispr*> crisprs = Array::get_crisprs(genome);
    start = time(start, "get crisprs");

    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    start = time(start, "get translations");

    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    start = time(start, "get fragments");
    
    start = time(start, genome_path.c_str());

    vector<Locus*> loci;
    for (Crispr* c : crisprs) loci.push_back(c);
    for (Fragment* f : fragments) loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) {return a->get_start() < b->get_start(); });

    write_loci(genome_path, loci);
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

string Fragment::to_string_debug()
{
    string amino_domain = Util::translate_genome(*reference_genome, genome_begin, genome_final, reference_translation->pos);
    string amino_gene = Util::translate_genome(*reference_genome, expanded_genome_begin, expanded_genome_final, reference_translation->pos);

    string dna_domain = reference_genome->substr(genome_begin, genome_final - genome_begin);
    string dna_gene = reference_genome->substr(expanded_genome_begin, expanded_genome_final - expanded_genome_begin);

    string amino_buffer;

    ui begin_discrepant = (genome_begin - expanded_genome_begin);
    ui final_discpreant = (expanded_genome_final - genome_final);

    dna_gene.insert(begin_discrepant, "-");
    amino_gene.insert(begin_discrepant / 3, "-");

    dna_gene.insert(dna_gene.size() - final_discpreant, "-");
    amino_gene.insert(amino_gene.size() - (final_discpreant / 3), "-");

    std::ostringstream out;
    out << fmt::format("{}\t{}\n", reference_translation->pos ? "+" : "-", reference_profile->gn);
    out << fmt::format("\t{}...{}\n", genome_begin, genome_final);
    out << fmt::format("\t{}...{}\n", expanded_genome_begin, expanded_genome_final);
    out << fmt::format("\t{}\n", amino_gene);
    out << fmt::format("\t{}\n", dna_gene);
    return out.str();
}

string Fragment::to_string_summary()
{
    string domain = reference_profile->gn;
    string strand = reference_translation->pos ? "+" : "-";
    string gn = domain_to_gn.contains(domain) ? domain_to_gn.at(domain) : domain;
    return fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", expanded_genome_begin, expanded_genome_final, strand, gn, genome_begin, genome_final, domain);
}

void debug()
{
    sage_interpreter("T:\\crispr\\groundtruth\\report_align.tsv", genome_dir);
    //string a = genome_dir + "\\GCA_000189575.1_ASM18957v1_genomic.fna";
    //prospect_genome(a);
    exit(0);
}

void close()
{
    out_debug.close();
    out_gene.close();
    out_domain.close();
}

int main()
{
    auto start_main = time();

    init();

    //debug();

    prospect_genome_dir(genome_dir);
    
    close();
    
    start_main = time(start_main, "main");
    return 0;                                                                                                           
}


