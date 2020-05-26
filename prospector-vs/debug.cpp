#include "debug.h"

//void Debug::map()
//{
    // ui query = 1283501;
    // ui q = genome_encoding[query];
    // for (ui i = 0; i < 1000; i++)
    // {
    //     ui pos = query + Prospector::k_start + Prospector::spacer_skip + i;
    //     ui diff = difference_cpu(genome_encoding[query], genome_encoding[pos]);

    //     printf("%s %d %d\n", genome.substr(pos, SIZE).c_str(), pos, diff);
    // }
//}


// visualize proximals
// for ( auto const& [query, proximal] : proximal_targets)
// {
//     // vector<ui> proximals = proximal_targets[q_i];
//     if (proximal.size() >= Prospector::repeats_min)
//     {
//         // ui query = queries[q_i];
//         fmt::print("{}:{}-{}\n", query, genome.substr(query, 16), genome.substr(query+16, 16));
//         for (ui target : proximal)
//         {
//             // ui qmap_index = q_i * Prospector::map_size_big + t_i;
//             // ui target = query + Prospector::k_start + Prospector::spacer_skip + t_i;
//             fmt::print("\t{}:{}-{}\n", target, genome.substr(target, 16), genome.substr(target+16, 16) );
//         }
//     }
//     fmt::print("\n");
// }



// void Debug::debug_clusters(const vector<vector<ull>>& clusters)
// {
//     for (vector<ull> cluster : clusters) 
//         fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
// }

vector<CasProfile> Debug::cas_filter(vector<CasProfile> profiles, string gn)
{
    return Util::filter(profiles, [&](CasProfile p) {return p.gn == gn; });
}

vector<Crispr>Debug::crispr_filter(vector<Crispr> crisprs, ui start, ui end)
{
    return Util::filter(crisprs, [&](Crispr c) { return c.start > start && c.end < end; });
}

string Debug::translation_test(const string& genome, ui genome_start, ui genome_final, bool pos, ui debug_aminos)
{

    genome_start -= debug_aminos * 3;
    genome_final += debug_aminos * 3;

    string translation = Util::translate_genome(genome, genome_start, genome_final, pos);

    auto len = translation.length();
    auto a = translation.substr(0, debug_aminos);
    auto b = translation.substr(debug_aminos, len - (debug_aminos * 2));
    auto c = translation.substr(len - debug_aminos);


    return fmt::format("\t\t\t{}--{}--{}\n", a, b, c);
}

void Debug::crispr_print(vector<Crispr> crisprs, const string& genome, ui start, ui end)
{
    auto filtered = crispr_filter(crisprs, start, end);
    int how_many = filtered.size();
    for (ull i = filtered.size() - how_many; i < filtered.size(); i++)
    {
        filtered[i].print(genome);
    }
    fmt::print("terminating after debug\n");
    exit(0);
}