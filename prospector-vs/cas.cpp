#include "cas.h"

static unordered_set<string> alternate_starts_pos { "GTG", "TTG" };
static unordered_set<string> alternate_starts_neg { "CAC", "CAA" };

static unordered_set<string> start_codons_pos{ "ATG" };
static unordered_set<string> start_codons_neg{ "CAT" };
static unordered_set<string> stop_codons_pos{ "TAA", "TAG", "TGA" };
static unordered_set<string> stop_codons_neg{ "TTA", "CTA", "TCA" };

// potential performance optimization: compare amino_family acids of codons in the amino_family rather than codons in the genome.
// should be faster to compare single chars rather than strings.

ll generic_expand(const string& genome, unordered_set<string>& sinks, ll begin, ll increment, ll increments)
{
    ll offset = 0;
    try
    {
        for (ll i = 0; i < increments; i++)
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

ll steps = 200;

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

vector<vector<ll>> cluster_index(const vector<ll>& indices)
{
    vector<vector<ll>> clusters;
    vector<ll> cluster;
    ll prev = indices[0];
    for (ll index : indices)
    {
        if (index - prev > Cas::max_inter_cluster_dist)
        {
            vector<ll> cluster_cp = cluster; 
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

bool good_clusters(const vector<vector<ll>>& clusters)
{
    for (vector<ll> cluster : clusters)
    {   
        if (cluster.size() > Cas::cluster_metric_min)
        {
            return true;
        }
    }
    return false;
}

ll demarc_start_clusters(const vector<vector<ll>>& clusters)
{
    for (const vector<ll>& cluster : clusters)
        if (cluster.size() > Cas::cluster_metric_min)
            return cluster[0];
    assert(false); return -1;
}

ll demarc_final_clusters(const vector<vector<ll>>& clusters)
{
    for (ll i = clusters.size()-1; i >= 0; i--)
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
    ll index_kmer_start = fragment->clust_begin;
    ll index_kmer_final = fragment->clust_final;

    // string protein = f.reference_translation->pure.substr(index_kmer_start, (index_kmer_end - index_kmer_start) + CasUtil::k);

    ll raw_pos_start = fragment->reference_translation->pure_mapping[index_kmer_start];
    ll raw_pos_end = fragment->reference_translation->pure_mapping[index_kmer_final];

    ll genome_begin;
    ll genome_final;
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

vector<Fragment*> Cas::cas(vector<CasProfile*>& profiles, vector<Translation*>& translations, string& genome)
{

    vector<Fragment*> fragments;

    #pragma omp parallel for
    for (signed long p = 0; p < profiles.size(); p++)
    {
        CasProfile* profile = profiles[p];
        for (Translation* t : translations)
        {
            vector<ll> index;

            for (ll i = 0; i < t->pure_kmerized_encoded.size(); i++)
            {
                bool contains = profile->hash_table.contains(t->pure_kmerized_encoded[i]);
                if (contains)
                    index.push_back(i);
            }

            if (index.size() == 0)
                continue;

            vector<vector<ll>> clusters = cluster_index(index);

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

            if (f->genome_begin < 0)
            {
                printf("break\n");
            }

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



vector<Translation*> Cas::get_triframe(const string& genome, ll genome_start, ll genome_final, bool pos)
{
    string domain = genome.substr(genome_start, genome_final - genome_start);
    //domain = pos ? domain : Util::reverse_complement(domain);

    if (!pos)
    {
        Util::reverse_complement(domain);
    }

    vector<Translation*> translations;
    for (ll frame = 0; frame < 3; frame++)
	{
        Translation* translation = new Translation;
        translation->pos = pos;
        translation->genome_start = pos ? genome_start + frame : genome_start;
        translation->genome_final = pos ? genome_final : genome_final - frame;

        translation->raw = Util::translate_domain(domain.substr(frame));

		translation->pure = "";

		ll stop_count = 0;
		ll index = 0;
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

vector <Translation*> Cas::get_sixframe(const string& genome, ll genome_start, ll genome_final)
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


    std::sort(crisprs.begin(), crisprs.end(), [](Crispr* a, Crispr* b) { return a->genome_start < b->genome_start; });
    const ll min_translation_size = 20;
    const ll permitted_translation_overlap = 100;
    ll prev_final = 0;
    for (Crispr* c : crisprs)
    {
        ll t_start = c->genome_start - Cas::upstream_size; // normal case
        ll t_final = c->genome_final + Cas::upstream_size; // normal case

        t_start = std::max(t_start, prev_final - permitted_translation_overlap); // previous translation may overlap

        // genome bounds
        t_start = std::max(t_start, (ll)0);
        t_final = std::min(t_final, (ll)genome.size() - 1);

        ll translation_size = t_final - t_start;

        if (translation_size < min_translation_size) continue;

        prev_final = t_final;

        // TODO: if the crispr is contained within the translation (it normally will be) then strip it out before computing the translation

        vector<Translation*> six_frame = Cas::get_sixframe(genome, t_start, t_final);
        for (Translation* t : six_frame) t->reference_crispr = c;
        translations.insert(translations.end(), six_frame.begin(), six_frame.end());
    }
    return translations;
}


bool Fragment::is_domain()
{
    bool is_domain = !CasProfileUtil::domain_table_contains(this->reference_profile->identifier);
    //fmt::print("is domain {} {}\n", this->reference_profile->identifier, is_domain);
    return is_domain;
}

bool Fragment::is_gene()
{
    return !this->is_domain();
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
    out << fmt::format("{}\t{}\n", reference_translation->pos ? "+" : "-", reference_profile->identifier);
    out << fmt::format("\t{}...{}\n", genome_begin, genome_final);
    out << fmt::format("\t{}...{}\n", expanded_genome_begin, expanded_genome_final);
    out << fmt::format("\t{}\n", amino_gene);
    out << fmt::format("\t{}\n", dna_gene);
    return out.str();
}

string Fragment::to_string_summary()
{
    string domain = reference_profile->identifier;
    string strand = reference_translation->pos ? "+" : "-";
    string identifier = CasProfileUtil::domain_table_contains(domain) ? CasProfileUtil::domain_table_fetch(domain) : domain;
    return fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", expanded_genome_begin, expanded_genome_final, strand, identifier, genome_begin, genome_final, domain);
}