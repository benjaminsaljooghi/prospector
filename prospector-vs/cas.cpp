#include "cas.h"
#include "debug.h"


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

bool* Cas::compute_target_map(const vector<const CasProfile*>& profiles, const vector<Translation*>& translations )
{
    auto start = time();

    ull num_translations = translations.size();
    ull num_cas = profiles.size();
    ull per_translation = (Cas::upstream_size - 1 + 3) / 3; // maximum size of a amino_family is CasUtil::upstream_size / 3 = 3333.33 (3334), but they are a bit smaller because of stop codons
    ull per_cas = per_translation * num_translations;


    ull target_map_size = num_cas * per_cas; 
    bool* target_map = (bool*) malloc(sizeof(bool) * target_map_size);

    #pragma omp parallel for
    for (signed long cas_i = 0; cas_i < num_cas; cas_i++)
    {  
        const CasProfile* cas_profile = profiles[cas_i];
        for (ull translation_i = 0; translation_i < num_translations; translation_i++)
        {
            for (ui i = 0; i < translations[translation_i]->pure_kmerized_encoded.size(); i++)
            {
                ui query = translations[translation_i]->pure_kmerized_encoded[i];
                string query_kmer = translations[translation_i]->pure_kmerized[i];

                ull target_map_index = (cas_i * per_cas) + (translation_i * per_translation) + i;

                // A
                //auto at_index = cas_profile.hash_table[query % cas_profile.N]; 
                //bool contains = at_index == query;

                // B
                bool contains = cas_profile->hash_table.contains(query);

                //fmt::print("{}:{}:{}:{}:{}:{}\n", target_map_index, cas_i, translation_i, i, contains, query_kmer, query);

                target_map[target_map_index] = contains;
            }
        }
    }

    start = time(start, "target map construction");

    return target_map;
}



static string start_codon_pos = "ATG";
static string start_codon_neg = "CAT";
static unordered_set<string> stop_codons_pos {"TAA", "TAG", "TGA"};
static unordered_set<string> stop_codons_neg {"TTA", "CTA", "TCA"};

// potential performance optimization: compare amino_family acids of codons in the amino_family rather than codons in the genome.
// should be faster to compare single chars rather than strings.
void fragment_expansion_pos(Fragment* fragment, const string& genome)
{
    fragment->expanded_genome_begin = fragment->genome_begin;
    fragment->expanded_genome_final = fragment->genome_final;
    for (ull i = 0; i < 500; i += 3)
    {
        ull new_begin = fragment->genome_begin - i;
        if (genome.substr(new_begin, 3) == start_codon_pos)
        {
            fragment->expanded_genome_begin = new_begin;
            break;
        }
    }

    for (ull i = 0; i < 500; i += 3)
    {
        ull new_final = fragment->genome_final + i;
        if (stop_codons_pos.contains(genome.substr(new_final - 3, 3)))
        {
            fragment->expanded_genome_final = new_final;
            break;
        }
    }
}

void fragment_expansion_neg(Fragment* fragment, const string& genome)
{
    fragment->expanded_genome_begin = fragment->genome_begin;
    fragment->expanded_genome_final = fragment->genome_final;
    for (ull i = 0; i < 500; i += 3)
    {
        ull new_begin = fragment->genome_begin - i;
        if (stop_codons_neg.contains(genome.substr(new_begin, 3)))
        {
            fragment->expanded_genome_begin = new_begin;
            break;
        }
    }

    for (ull i = 0; i < 500; i += 3)
    {
        ull new_final = fragment->genome_final + i;
        if (start_codon_neg == genome.substr(new_final - 3, 3))
        {
            fragment->expanded_genome_final = new_final;
            break;
        }
    }

}


vector<Fragment*> Cas::cas(const vector<const CasProfile*>& profiles, const vector<Translation*>& translations, const string& genome)
{

    fmt::print("\tidentifying cas genes in {} translations...\n", translations.size());
    auto start = time();  

    bool* target_map = compute_target_map(profiles, translations);

    ull num_translations = translations.size();
    ull num_cas = profiles.size();
    ull per_translation = 3334; // maximum size of a amino_family is CasUtil::upstream_size / 3 = 3333.33 (3334), but they are a bit smaller because of stop codons
    ull per_cas = per_translation * num_translations;

    vector<Fragment*> fragments;
    for (ull cas_i = 0; cas_i < num_cas; cas_i++)
    {
        for (ull translation_i = 0; translation_i < num_translations; translation_i++)
        {
            vector<ull> index;
            for (ull i = 0; i < translations[translation_i]->pure_kmerized_encoded.size(); i++)
            {
                ull target_map_index = (cas_i * per_cas) + (translation_i * per_translation) + i;
                if (target_map[target_map_index])
                {
                    index.push_back(i);
                }
            }


            if (index.size() == 0)
                continue;

            vector<vector<ull>> clusters = cluster_index(index);


            // fmt::print("{}:{}\n", cas_i, translation_i);
            // debug_clusters(clusters);

            if (!good_clusters(clusters))
                continue;

            Fragment* f = new Fragment;
            f->reference_crispr = translations[translation_i]->reference_crispr;
            f->reference_translation = translations[translation_i];
            f->reference_profile = profiles[cas_i];
            f->clusters = clusters;

            compute_demarc(f);

            if (f->clust_final - f->clust_begin <= 15)
            {
                continue;
            }

            compute_details(f, genome);

            if (f->reference_translation->pos)
            {
                fragment_expansion_pos(f, genome);
            }
            else
            {
                fragment_expansion_neg(f, genome);
            }

            fragments.push_back(f);
           
        }
    }
    start = time(start, "target map parse");


    //for (Fragment& a : fragments) compute_demarc(a);
    //genome_start = time(genome_start, "f demarcs");

    //for (Fragment& a : fragments) compute_details(a, genome);
    //genome_start = time(genome_start, "f details");

    return fragments;
}

Gene* gene_from_fragments(vector<Fragment*>& fragments)
{
    sort(fragments.begin(), fragments.end(), [](const Fragment* a, const Fragment* b) { return a->genome_begin < b->genome_begin;  });

    Gene* g = new Gene;
    g->gn = fragments[0]->reference_profile->gn;
    g->fragments = fragments;
    return g;
}

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

vector<Gene*> genes_from_fragments(vector<Fragment*>& fragments)
{
    map<string, vector<Fragment*>> gene_fragments;
    for (Fragment* a : fragments)
        gene_fragments[a->reference_profile->gn].push_back(a);

    vector<Gene*> genes;
    for (auto gene : gene_fragments)
    {
        genes.push_back(gene_from_fragments(gene.second));
    }

    //genes = best_genes(genes);
    sort(genes.begin(), genes.end(), [](Gene* a, Gene* b) { return a->fragments[0]->genome_begin < b->fragments[0]->genome_begin; } );

    return genes;
}

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

string crispr_type(vector<Gene*>& genes)
{
    for (Gene* gene : genes)
    {
        if (type_lookup.contains(gene->gn))
        {
            return type_lookup.at(gene->gn);
        }
    }
    return "?";
}

map<string, vector<Gene*>> Cas::assemble_genes(const vector<Crispr*>& crisprs, const vector<Fragment*>& fragments)
{
    //map<string, Crispr> crispr_map;

    map<string, vector<Fragment*>> crispr_fragments;
    map<string, vector<Gene*>> crispr_genes;


    for (const Crispr* c : crisprs)
    {
        string c_string = c->identifier_string();
        for (Fragment* f : fragments)
        {
            if (f->reference_crispr->start == c->start && f->reference_crispr->k == c->k)
            {
                crispr_fragments[c_string].push_back(f); 
            }
        }
    }

    for (auto [c_string, fragments] : crispr_fragments)
    {
        crispr_genes[c_string] = genes_from_fragments(fragments);
    }

    return crispr_genes;
}

//void print_gene_summary(Gene& gene)
//{
//    fmt::print("\t{}\n", gene.reference_profile->gn);
//    for (const Fragment& f : gene.fragments)
//    {
//        fmt::print("\t\t{}...{}\n", f.details->genome_start, f.details->genome_final);
//    }
//    fmt::print("\n");
//}

void Cas::print_fragment_debug(const Fragment* f, const string& genome)
{

    string amino_family = Util::translate_genome(genome, f->genome_begin, f->genome_final, f->reference_translation->pos);
    string amino_cds = Util::translate_genome(genome, f->expanded_genome_begin, f->expanded_genome_final, f->reference_translation->pos);

    string dna_family = genome.substr(f->genome_begin, f->genome_final - f->genome_begin);
    string dna_cds = genome.substr(f->expanded_genome_begin, f->expanded_genome_final - f->expanded_genome_begin);

    fmt::print("\t{}\n", f->reference_translation->reference_crispr->identifier_string());
    //fmt::print("\t{}, \n", f->reference_crispr->end - f->reference_translation->genome_start, f->refe);

    fmt::print("\t{}...{}\n", f->genome_begin, f->genome_final);
    fmt::print("\t{}...{}\n", f->expanded_genome_begin, f->expanded_genome_final);

    fmt::print("\t{}\n", amino_family);
    fmt::print("\t{}\n", amino_cds);

    fmt::print("\t{}\n", dna_family);
    fmt::print("\t{}\n", dna_cds);
}

void print_gene_debug(Gene* gene, const string& genome)
{
    //ui size = gene->size();
    fmt::print("\t{}:{}\n", gene->gn, gene->fragments[0]->reference_translation->pos);
    for (const Fragment* fragment : gene->fragments)
    {
        Cas::print_fragment_debug(fragment, genome);
    }

    fmt::print("\n");
}


void Cas::print_all(const vector<Crispr*>& crisprs, const map<string, vector<Gene*>>& crispr_genes, const string& genome)
{
    for (const Crispr* c : crisprs)
    {
        c->print(genome);
        string c_string = c->identifier_string();
        if (!crispr_genes.contains(c_string))
        {
            continue;
        }

        vector<Gene*> genes = crispr_genes.at(c_string);
        string _type = crispr_type(genes);
        string _class = class_lookup.at(_type);

        for (Gene* gene : genes)
        {
            print_gene_debug(gene, genome);
            //print_gene_summary(gene);
        }

        fmt::print("{}:{}:{}\n\n", _class, _type, c_string);
    }
}


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


vector<Translation*> Cas::crispr_proximal_translations(const string& genome, const vector<Crispr*>& crisprs)
{
    auto start = time();
    vector<Translation*> translations;
    for (const Crispr* c : crisprs)
    {
        ull genome_start = c->start - Cas::upstream_size; genome_start = genome_start < c->start ? genome_start : 0;
        vector<Translation*> up = get_sixframe(genome, genome_start, c->start);
        
        ull genome_end = c->end + Cas::upstream_size; genome_end = genome_end > c->end ? genome_end : genome.size() - 1;
        vector<Translation*> down = get_sixframe(genome, c->end, genome_end);

        for (Translation* t : up)
        {
            t->reference_crispr = c;
            translations.push_back(t);
        }

        for (Translation* t : down)
        {
            t->reference_crispr = c;
            translations.push_back(t);
        }
    }
    time(start, "get translations");
    return translations;
}
