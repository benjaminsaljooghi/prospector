#include "debug.h"

void Debug::visualize_map(string& genome_path)
{
    string genome = Util::load_genome(genome_path);
    Prospector::Encoding encoding = Prospector::get_genome_encoding(genome.c_str(), genome.size());

    ui query = 2110970;
    auto k = 32; 
    for (ui i = 0; i < 50000; i++)
    {
        ui target = query + i;
        
        auto diff = Util::difference_cpu(encoding.h[query], encoding.h[target]);
        auto query_str = genome.substr(query, k);
        auto target_str = genome.substr(target, k);
        auto mutant = Array::mutant(genome.c_str(), encoding.h, k, query, target, Prospector::repeat_tolerance_ratio_sensitive);

        if (mutant)
            fmt::print("{} -> {} {} {} {} {}\n", query_str, target, target_str, target + k, diff, mutant);
    }
    exit(0);
}

void Debug::visualize_proximals(map<ull, vector<ull>> proximal_targets, string genome)
{
    for (auto const& [query, proximal] : proximal_targets)
    {
        // vector<ui> proximals = proximal_targets[q_i];
        //if (proximal.size() >= Prospector::repeats_min)
        {
            // ui query = queries[q_i];
            fmt::print("{}:{}-{}\n", query, genome.substr(query, 16), genome.substr(query + 16, 16));
            for (ull target : proximal)
            {
                // ui qmap_index = q_i * Prospector::map_size_big + t_i;
                // ui target = query + Prospector::k_start + Prospector::spacer_skip + t_i;
                fmt::print("\t{}:{}-{}\n", target, genome.substr(target, 16), genome.substr(target + 16, 16));
            }
        }
        fmt::print("\n");
    }
    exit(0);
}

// void Debug::debug_clusters(const vector<vector<ull>>& clusters)
// {
//     for (vector<ull> cluster : clusters) 
//         fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
// }

vector<CasProfile*> Debug::cas_filter(vector<CasProfile*> profiles, string identifier)
{
    return Util::filter(profiles, [&](CasProfile* p) {return p->identifier == identifier; });
}

vector<Crispr*> Debug::crispr_filter(vector<Crispr*> crisprs, ull start, ull end)
{
    return Util::filter(crisprs, [&](Crispr* c) { return c->genome_start > start && c->genome_final < end; });
}

void Debug::genome_substr(const string& genome_path, ull genome_start, ull genome_final)
{
    fmt::print("{}\n", Util::load_genome(genome_path).substr(genome_start, genome_final - genome_start));
    exit(0); 
}

string Debug::translation_test(const string& genome, ull genome_start, ull genome_final, bool pos, ull debug_aminos)
{
    genome_start -= debug_aminos * 3;
    genome_final += debug_aminos * 3;

    string translation = Util::translate_genome(genome, genome_start, genome_final, pos);

    auto len = translation.length();
    auto a = translation.substr(0, debug_aminos);
    auto b = translation.substr(debug_aminos, len - (debug_aminos * 2));
    auto c = translation.substr(len - debug_aminos);

    return fmt::format("{}--{}--{}", a, b, c);
}

void Debug::translation_print(const string& genome_path, ull genome_start, ull genome_final, bool pos, ull debug_aminos)
{
    string genome = Util::load_genome(genome_path);
    string translation = Debug::translation_test(genome, genome_start, genome_final, pos, debug_aminos);
    fmt::print("debug: {}..{} {}\n", genome_start, genome_final, translation);
    exit(0);
}

void Debug::triframe_print(const string& genome, ull genome_start, ull genome_final, bool pos)
{
    auto triframe = Cas::get_triframe(genome, genome_start, genome_final, pos);
    for (auto translation : triframe)
    {
        fmt::print("{}\n", translation->pure);
    }
    //return triframe;
}

void Debug::cas_detect(const string& genome_path, ull genome_start, ull genome_final, bool pos, CasProfile* profile)
{
    string genome = Util::load_genome(genome_path);
    Translation* translation = Cas::get_triframe(genome, genome_start, genome_final, pos)[0];
    
    fmt::print("translation raw: {}\n", translation->raw);
    
    fmt::print("containment info:\n");
    for (auto kmer : translation->pure_kmerized)
    {
        auto enco = Util::encode_amino_kmer(kmer);
        bool contains = profile->hash_table.contains(enco);
        fmt::print("{} : {}\n", kmer, contains);
    }

    vector<CasProfile*> profiles;
    vector<Translation*> translations;
    
    profiles.push_back(profile);
    translations.push_back(translation);
    
    vector<Fragment*> fragments = Cas::cas(profiles, translations, genome);
    fmt::print("fragment info:\n");
    for (int i = 0; i < fragments.size(); i++)
        fmt::print("{}\n", fragments[i]->to_string_debug());
}

void Debug::crispr_print(vector<Crispr*> crisprs, const string& genome, ull start, ull end)
{
    auto filtered = Debug::crispr_filter(crisprs, start, end);
    Util::sort(filtered, CrisprUtil::heuristic_greater);
    for (ull i = 0; i < filtered.size(); i++)
        fmt::print("{}\n", filtered[i]->to_string_debug());
    exit(0);
}




void Debug::cartograph_interpreter(std::filesystem::path path, std::filesystem::path genome_dir)
{
    std::ifstream in(path);

	if (!in.good())
		throw runtime_error("input not good!");
        
    string line;
    string genome_accession = "";
    string genome = "";

    string interpretation_path = "cartograph_interpretation.txt";
    std::ofstream interpretation(interpretation_path.c_str());

    auto gen_debug_str = [&genome](string domains, string signals, ull begin, ull final, string strand) {
        auto a = genome.substr(begin, final - begin);
        auto b = domains == "CRISPR" ? "" : Debug::translation_test(genome, begin, final, strand == "+", 0);

        std::ostringstream stream;
        stream << fmt::format("\t{}\t{}\t{}..{}\t{}\n", domains, signals, begin, final, strand);
        stream << "\t" << a << "\n";
        stream << "\t" << b << "\n";
        return stream.str();
    };

    auto gen_ground_str = [&](std::vector<string> split) {
        ull begin = std::stoull(split[1]);
        ull final = std::stoull(split[2]);
        string strand = split[3];
        string domains = split[4];
        string signals = split[5];
        return gen_debug_str(domains, signals, begin, final, strand);
    };

    auto gen_target_str = [&](std::vector<string> split) {
        ull begin = std::stoull(split[6]);
        ull final = std::stoull(split[7]);
        string strand = split[8];
        string domains = split[9];
        string signals = split[10];
        return gen_debug_str(domains, signals, begin, final, strand);
    };

    while (std::getline(in, line))
    {
        if (line == "" || line[0] == '-')
            continue;

        auto split = Util::parse(line, "\t");
        string alignment_type = split[0];        

        if (alignment_type == "===")
        {
            genome_accession = split[1];
            interpretation << fmt::format("=== {}\n", genome_accession);
            
            for (const auto& entry : filesystem::directory_iterator(genome_dir))
            {
                string genome_path = entry.path().string();
                if (genome_path.find(genome_accession) != string::npos)
                {
                    genome = Util::load_genome(genome_path);
                    break;
                }
            }

        }

        if (alignment_type == "<")
        {
            interpretation << "ground only:\n";
            interpretation << gen_ground_str(split);
        }
        
        if (alignment_type == ">")
        {
            interpretation << "target only:\n";
            interpretation << gen_target_str(split);
        }

        if (alignment_type == "!")
        {

        }

        if (alignment_type == "$")
        {

        }

        interpretation << "\n\n";
    }


    interpretation.close();

    fmt::print("wrote {}\n", interpretation_path);
}