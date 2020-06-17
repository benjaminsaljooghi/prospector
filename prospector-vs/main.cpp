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
std::ofstream out_debug("out_debug.txt");
std::ofstream out_results("out_results.txt");

std::vector<std::string> parse(std::string str, std::string delim)
{
    //std::string s = "scott>=tiger>=mushroom";
    //std::string delimiter = ">=";

    std::vector < std::string> tokens;

    size_t pos = 0;
    std::string token;
    while ((pos = str.find(delim)) != std::string::npos) {
        token = str.substr(0, pos);
        //std::cout << token << std::endl;
        str.erase(0, pos + delim.length());
        tokens.push_back(token);
    }
    return tokens;
    //std::cout << s << std::endl;
}

void sage_interpreter(string path, string genome_dir)
{
    std::ifstream in(path);

    string line;
    
    while (std::getline(in, line))
    {
        if (line[0] == '-')
        {
            break;
        }
    }
    
    
    string genome_accession = "";
    bool invalidate = false;
    string genome = "";

    std::ofstream sage_interpretation("sage_interpretation.txt");

    while (std::getline(in, line))
    {
        if (line[0] == '-')
        {
            continue;
        }

        auto split = parse(line, "\t");

        if (split[0] == "===")
        {
            bool invalidate = true;
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
        


        auto get_debug_str = [&genome](string kind, ui begin, ui final, string strand) {
            auto a = genome.substr(begin, final - begin);

            auto b = kind == "CRISPR" ? "" : Debug::translation_test(genome, begin, final, strand == "+", 0);
        
            std::ostringstream stream;

            stream << fmt::format("\t{}\t{}..{}\t{}\n", kind, begin, final, strand);
            stream << "\t" << a << "\n";
            stream << "\t" << b << "\n";

            return stream.str();
        };

        //std::ostringstream stream;

        

        if (split[2] != "")
        {
            ui g_begin = stoi(split[2]);
            ui g_final = stoi(split[3]);
            string g_strand = split[4];
            auto str = get_debug_str(split[0], g_begin, g_final, g_strand);
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
            auto str = get_debug_str(split[0], t_begin, t_final, t_strand);
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

void prospect_genome(string genome_path)
{
    fmt::print("\n\n\n");

    auto start = time();

    string genome = Util::load_genome(genome_path);

    vector<Crispr*> crisprs = Array::get_crisprs(genome);

    vector<Translation*> translations = Cas::crispr_proximal_translations(genome, crisprs);
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    
    start = time(start, genome_path.c_str());

    vector<Locus*> loci;
    for (Crispr* c : crisprs) loci.push_back(c);
    for (Fragment* f : fragments) loci.push_back(f);

    std::sort(loci.begin(), loci.end(), [](Locus* a, Locus* b) {return a->get_start() < b->get_start(); });

    out_results << get_header(genome_path);
    for (Locus* l : loci) out_results << l->to_string_summary();
    
    out_debug << get_header(genome_path);
    for (Locus* l : loci) out_debug << l->to_string_debug();
    
}

void prospect_genome_dir(string genome_dir)
{
    ui i = 0;
    for (const auto& entry : filesystem::directory_iterator(genome_dir))
    {
        if (i++ > 3) break;
        string genome_path = entry.path().string();
        prospect_genome(genome_path);
    }
    //return results;
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
    string genome_dir = "T:\\crispr\\supp\\genomes";
    Prospector::device_init();
    profiles = CasProfileUtil::deserialize("T:\\crispr\\cas\\serial");
    

    // ----------- debug --------
    //string debug_path = "T:\\crispr\\supp\\genomes\\GCA_000730285.1_ASM73028v1_genomic.fna";
    
    //Debug::visualize_map(debug_path);
    //prospect_genome(debug_path, results);

    sage_interpreter("T:\\crispr\\supp\\stats.tsv", genome_dir);

    // ----------- prospect --------------
    //prospect_genome_dir(genome_dir);
    

    // --------- close -----------
    out_debug.close();
    out_results.close();
    start_main = time(start_main, "main"); return 0;                                                                                                           
}


