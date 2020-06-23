#include "cas.h"

static unordered_set<string> alternate_starts_pos { "GTG", "TTG" };
static unordered_set<string> alternate_starts_neg { "CAC", "CAA" };

static unordered_set<string> start_codons_pos{ "ATG" };
static unordered_set<string> start_codons_neg{ "CAT" };
static unordered_set<string> stop_codons_pos{ "TAA", "TAG", "TGA" };
static unordered_set<string> stop_codons_neg{ "TTA", "CTA", "TCA" };

// potential performance optimization: compare amino_family acids of codons in the amino_family rather than codons in the genome.
// should be faster to compare single chars rather than strings.

ull generic_expand(const string& genome, unordered_set<string>& sinks, ull begin, ull increment, ull increments)
{
    ull offset = 0;
    try
    {
        for (ull i = 0; i < increments; i++)
        {
            if (sinks.contains(genome.substr(begin + offset, 3)))
            {
                return begin + offset;
            }
            offset += increment;
        }
    }
    catch (exception& e)
    {
        return begin;
    }
    return begin;
}

ull steps = 200;

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

    // string protein = f.reference_translation->pure.substr(index_kmer_start, (index_kmer_end - index_kmer_start) + CasUtil::k);

    ull raw_pos_start = fragment->reference_translation->pure_mapping[index_kmer_start];
    ull raw_pos_end = fragment->reference_translation->pure_mapping[index_kmer_final];

    ull genome_begin;
    ull genome_final;
    if (fragment->reference_translation->pos)
    {
        genome_begin = fragment->reference_translation->genome_start + (raw_pos_start * 3);
        genome_final = fragment->reference_translation->genome_start + ((raw_pos_end + CasProfileUtil::k) * 3) + 3;
    }
    else
    {
        genome_final = fragment->reference_translation->genome_final - (raw_pos_start * 3);
        genome_begin = fragment->reference_translation->genome_final - ( ((raw_pos_end + CasProfileUtil::k) * 3) + 3 );
    }

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

vector<Fragment*> Cas::cas(vector<CasProfile*>& profiles, vector<Translation*>& translations, string& genome)
{

    vector<Fragment*> fragments;

    #pragma omp parallel for
    for (signed long p = 0; p < profiles.size(); p++)
    {
        CasProfile* profile = profiles[p];
        for (Translation* t : translations)
        {
            vector<ull> index;

            for (ull i = 0; i < t->pure_kmerized_encoded.size(); i++)
            {
                bool contains = profile->hash_table.contains(t->pure_kmerized_encoded[i]);
                if (contains)
                    index.push_back(i);
            }

            if (index.size() == 0)
                continue;

            vector<vector<ull>> clusters = cluster_index(index);

            if (!good_clusters(clusters))
                continue;

            Fragment* f = new Fragment;
            f->reference_genome = &genome;
            f->reference_crispr = t->reference_crispr;
            f->reference_translation = t;
            f->reference_profile = profile;
            f->clusters = clusters;

            compute_demarc(f);

            if (f->clust_final - f->clust_begin <= 15)
                continue;

            compute_details(f, genome);

            auto expansion = f->reference_translation->pos ? fragment_expansion_pos : fragment_expansion_neg;
            expansion(f, genome);

            #pragma omp critical
            {
                fragments.push_back(f);
            }
        }
    }


    return fragments;
}

//Gene* gene_from_fragments(vector<Fragment*>& fragments)
//{
//    sort(fragments.begin(), fragments.end(), [](const Fragment* a, const Fragment* b) { return a->genome_begin < b->genome_begin;  });
//
//    Gene* g = new Gene;
//    g->gn = fragments[0]->reference_profile->gn;
//    g->fragments = fragments;
//    return g;
//}

//vector<Gene> best_genes(vector<Gene> genes)
//{
//    map<string, Gene> types;
//    for (const Gene& gene : genes)
//    {
//        auto gn = gene.fragments[0].reference_profile->gn;
//        if (types.contains(gn) && gene.size() <= types.at(gn).size())
//            continue;
//        types[gn] = gene;
//    }
//    vector<Gene> __genes;
//    for (auto pairing : types)
//        __genes.push_back(pairing.second);
//        
//    return __genes;
//}

//vector<Gene*> genes_from_fragments(vector<Fragment*>& fragments)
//{
//    map<string, vector<Fragment*>> gene_fragments;
//    for (Fragment* a : fragments)
//        gene_fragments[a->reference_profile->gn].push_back(a);
//
//    vector<Gene*> genes;
//    for (auto gene : gene_fragments)
//    {
//        genes.push_back(gene_from_fragments(gene.second));
//    }
//
//    //genes = best_genes(genes);
//    sort(genes.begin(), genes.end(), [](Gene* a, Gene* b) { return a->fragments[0]->genome_begin < b->fragments[0]->genome_begin; } );
//
//    return genes;
//}


//map<string, vector<Gene*>> Cas::assemble_genes(const vector<Crispr*>& crisprs, const vector<Fragment*>& fragments)
//{
//    //map<string, Crispr> crispr_map;
//
//    map<string, vector<Fragment*>> crispr_fragments;
//    map<string, vector<Gene*>> crispr_genes;
//
//
//    for (const Crispr* c : crisprs)
//    {
//        string c_string = c->identifier_string();
//        for (Fragment* f : fragments)
//        {
//            if (f->reference_crispr->start == c->start && f->reference_crispr->k == c->k)
//            {
//                crispr_fragments[c_string].push_back(f); 
//            }
//        }
//    }
//
//    for (auto [c_string, fragments] : crispr_fragments)
//    {
//        crispr_genes[c_string] = genes_from_fragments(fragments);
//    }
//
//    return crispr_genes;
//}



vector<Translation*> Cas::get_triframe(const string& genome, ull genome_start, ull genome_final, bool pos)
{
    string domain = genome.substr(genome_start, genome_final - genome_start);
    //domain = pos ? domain : Util::reverse_complement(domain);

    if (!pos)
    {
        Util::reverse_complement(domain);
    }

    vector<Translation*> translations;
    for (ull frame = 0; frame < 3; frame++)
	{
        Translation* translation = new Translation;
        translation->pos = pos;
        translation->genome_start = pos ? genome_start + frame : genome_start;
        translation->genome_final = pos ? genome_final : genome_final - frame;

        translation->raw = Util::translate_domain(domain.substr(frame));

		translation->pure = "";

		ull stop_count = 0;
		ull index = 0;
		for (char elem : translation->raw )
		{
			if (elem == Util::stop_c)
			{
				stop_count++;
				continue;
			}
			translation->pure += elem;
			translation->pure_mapping.push_back(index + stop_count);
			index++;
		}
        translation->pure_kmerized = Util::kmerize(translation->pure, CasProfileUtil::k);
        translation->pure_kmerized_encoded = Util::encode_amino_kmers(translation->pure_kmerized, CasProfileUtil::k);
        translations.push_back(translation);
    }
    return translations;
}

vector <Translation*> Cas::get_sixframe(const string& genome, ull genome_start, ull genome_final)
{
    vector<Translation*> sixframe;

    for (Translation* translation : Cas::get_triframe(genome, genome_start, genome_final, true))
        sixframe.push_back(translation);

    for (Translation* translation : Cas::get_triframe(genome, genome_start, genome_final, false))
        sixframe.push_back(translation);

    return sixframe;
}


vector<Translation*> Cas::crispr_proximal_translations(const string& genome, vector<Crispr*>& crisprs)
{
    vector<Translation*> translations;

    ull max = 0;

    std::sort(crisprs.begin(), crisprs.end(), [](Crispr* a, Crispr* b) { return a->start < b->start; });

    for (signed long c_index = 0; c_index < crisprs.size(); c_index++)
    {
        Crispr* c = crisprs[c_index];

        ull l_begin = c->start - (ull) Cas::upstream_size; 
        ull l_final = c->start;

        ull r_begin = c->end;
        ull r_final = c->end + (ull) Cas::upstream_size;

        if (l_begin > c->start) // overflow
            l_begin = 0;

        if (r_final >= genome.size())
            r_final = genome.size() - 1;

        l_begin = std::max(max, l_begin);
        r_begin = std::max(max, r_begin);

        bool can_l = l_begin < l_final;
        bool can_r = r_begin < genome.size() && r_begin < r_final;

        if (can_r)
        {
            max = std::max(max, r_final);
        }
        
        vector<Translation*> l;
        vector<Translation*> r;

        if (can_l)
            l = get_sixframe(genome, l_begin, l_final);

        if (can_r)
            r = get_sixframe(genome, r_begin, r_final);

        vector<Translation*> local;
        local.insert(local.end(), l.begin(), l.end());
        local.insert(local.end(), r.begin(), r.end());
        for (Translation* t : local) t->reference_crispr = c;

        translations.insert(translations.end(), local.begin(), local.end());
        
    }
    return translations;
}


