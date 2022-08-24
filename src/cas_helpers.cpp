#include "cas_helpers.h"

ull generic_expand(const string& genome, unordered_set<string>& sinks, ull begin, ull increment, ull increments)
{
    ull offset = 0;
    try
    {
        for (ull i = 0; i < increments; i++)
        {
            if (sinks.find(genome.substr(begin + offset, 3)) != sinks.end())
            {
                return begin + offset;
            }
            offset += increment;
        }
    }
    catch (exception&)
    {
        return begin;
    }
    return begin;
}


void fragment_expansion_pos(Fragment* fragment, const string& genome)
{
    fragment->expanded_genome_begin = generic_expand(genome, start_codons_pos, fragment->genome_begin, -3, steps);
    fragment->expanded_genome_final = generic_expand(genome, stop_codons_pos, fragment->genome_final - 3, 3, steps) + 3;

    if (fragment->expanded_genome_begin == fragment->genome_begin)
    {
        fragment->expanded_genome_begin = generic_expand(genome, alternate_starts_pos, fragment->genome_begin, -3, steps);
    }
}

void fragment_expansion_neg(Fragment* fragment, const string& genome)
{
    fragment->expanded_genome_begin = generic_expand(genome, stop_codons_neg, fragment->genome_begin, -3, steps);
    fragment->expanded_genome_final = generic_expand(genome, start_codons_neg, fragment->genome_final - 3, 3, steps) + 3;

    if (fragment->expanded_genome_final == fragment->genome_final)
    {
        fragment->expanded_genome_final = generic_expand(genome, alternate_starts_neg, fragment->genome_final - 3, 3, steps) + 3;
    }

}

vector<vector<ull>> cluster_index(const vector<ull>& indices)
{
    vector<vector<ull>> clusters;
    vector<ull> cluster;
    ull prev = indices[0];
    for (ull index : indices)
    {
        if (index - prev > Cas::max_inter_cluster_dist)
        {
            vector<ull> cluster_cp = cluster;
            if (cluster_cp.size() > 1)
            {
                clusters.push_back(cluster_cp);
            }
            cluster.clear();
        }
        cluster.push_back(index);
        prev = index;
    }
    clusters.push_back(cluster);
    return clusters;
}

bool good_clusters(const vector<vector<ull>>& clusters)
{
    for (vector<ull> cluster : clusters)
    {
        if (cluster.size() > Cas::cluster_metric_min)
        {
            return true;
        }
    }
    return false;
}

ull demarc_start_clusters(const vector<vector<ull>>& clusters)
{
    for (const vector<ull>& cluster : clusters)
        if (cluster.size() > Cas::cluster_metric_min)
            return cluster[0];
    assert(false); return -1;
}

ull demarc_final_clusters(const vector<vector<ull>>& clusters)
{
    for (ull i = clusters.size()-1; i >= 0; i--)
        if (clusters[i].size() > Cas::cluster_metric_min)
            return clusters[i][clusters[i].size()-1];
    assert(false); return -1;
}

void compute_demarc(Fragment* frag)
{
    frag->clust_begin = demarc_start_clusters(frag->clusters);
    frag->clust_final = demarc_final_clusters(frag->clusters);
}

void compute_details(Fragment* fragment, const string& genome)
{
    ull index_kmer_start = fragment->clust_begin;
    ull index_kmer_final = fragment->clust_final;

    fragment->protein_sequence = fragment->reference_translation->pure.substr(index_kmer_start, (index_kmer_final - index_kmer_start) + CasProfileUtil::k);

    ull raw_pos_start = fragment->reference_translation->pure_mapping[index_kmer_start];
    ull raw_pos_end = fragment->reference_translation->pure_mapping[index_kmer_final];

    ull genome_begin;
    ull genome_final;
    if (fragment->reference_translation->pos)
    {
        genome_begin = fragment->reference_translation->genome_start + (raw_pos_start * 3); // inclusive begin
        genome_final = fragment->reference_translation->genome_start + ((raw_pos_end + CasProfileUtil::k) * 3); // exclusive final
    }
    else
    {
        genome_final = fragment->reference_translation->genome_final - (raw_pos_start * 3); // exclusive final
        genome_begin = fragment->reference_translation->genome_final - ( ((raw_pos_end + CasProfileUtil::k) * 3)); // inclusive begin
    }

    genome_final += 1; // to be exclusive.

    //string genome_translation = Debug::translation_test(genome, genome_start, genome_final, f.reference_translation->pos, 6);
    //string amino_family = Util::translate_genome(genome, genome_start, genome_final, f.reference_translation->pos);


    fragment->genome_begin = genome_begin;
    fragment->genome_final = genome_final;
    //f->genome_translation = genome_translation;
    //details->quality = genome_final - genome_start;

    //f.details = details;
}

bool fragment_equivalent(const Fragment* a, const Fragment* b)
{
    return a->clust_begin == b->clust_begin && a->clust_final == b->clust_final;
}

bool fragment_contains(const Fragment* a, const Fragment* b)
{
    return a->clust_begin > b->clust_begin && a->clust_final < b->clust_final;
}

