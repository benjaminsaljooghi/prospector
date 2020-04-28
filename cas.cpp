#include "cas.h"

size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2)
{
  const size_t m(s1.size());
  const size_t n(s2.size());
 
  if( m==0 ) return n;
  if( n==0 ) return m;
 
  size_t *costs = new size_t[n + 1];
 
  for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
  size_t i = 0;
  for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
  {
    costs[0] = i+1;
    size_t corner = i;
 
    size_t j = 0;
    for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
    {
      size_t upper = costs[j+1];
      if( *it1 == *it2 )
      {
		  costs[j+1] = corner;
	  }
      else
	  {
		size_t t(upper<corner?upper:corner);
        costs[j+1] = (costs[j]<t?costs[j]:t)+1;
	  }
 
      corner = upper;
    }
  }
 
  size_t result = costs[n];
  delete [] costs;
 
  return result;
}


const map <string, string> codon_table = {
    {"TTT", "F"},
    {"TTC", "F"},
    {"TTA", "L"},
    {"TTG", "L"},
    {"CTT", "L"},
    {"CTC", "L"},
    {"CTA", "L"},
    {"CTG", "L"},
    {"ATT", "I"},
    {"ATC", "I"},
    {"ATA", "I"},
    {"ATG", "M"},
    {"GTT", "V"},
    {"GTC", "V"},
    {"GTA", "V"},
    {"GTG", "V"},
    {"TCT", "S"},
    {"TCC", "S"},
    {"TCA", "S"},
    {"TCG", "S"},
    {"CCT", "P"},
    {"CCC", "P"},
    {"CCA", "P"},
    {"CCG", "P"},
    {"ACT", "T"},
    {"ACC", "T"},
    {"ACA", "T"},
    {"ACG", "T"},
    {"GCT", "A"},
    {"GCC", "A"},
    {"GCA", "A"},
    {"GCG", "A"},
    {"TAT", "Y"},
    {"TAC", "Y"},
    {"TAA", CasUtil::stop},
    {"TAG", CasUtil::stop},
    {"CAT", "H"},
    {"CAC", "H"},
    {"CAA", "Q"},
    {"CAG", "Q"},
    {"AAT", "N"},
    {"AAC", "N"},
    {"AAA", "K"},
    {"AAG", "K"},
    {"GAT", "D"},
    {"GAC", "D"},
    {"GAA", "E"},
    {"GAG", "E"},
    {"TGT", "C"},
    {"TGC", "C"},
    {"TGA", CasUtil::stop},
    {"TGG", "W"},
    {"CGT", "R"},
    {"CGC", "R"},
    {"CGA", "R"},
    {"CGG", "R"},
    {"AGT", "S"},
    {"AGC", "S"},
    {"AGA", "R"},
    {"AGG", "R"},
    {"GGT", "G"},
    {"GGC", "G"},
    {"GGA", "G"},
    {"GGG", "G"}
};

map<char, ui> amino_encoding {
    {'F', 0},
    {'L', 1},
    {'I', 2},
    {'M', 3},
    {'V', 4},
    {'S', 5},
    {'P', 6},
    {'T', 7},
    {'A', 8},
    {'Y', 9},
    {'H', 10},
    {'Q', 11},
    {'N', 12},
    {'K', 13},
    {'D', 14},
    {'E', 15},
    {'C', 16},
    {'W', 17},
    {'R', 18},
    {'G', 19},
};

string translate_domain(const string& domain)
{
    ull codon_size = 3;
    string raw = "";
    for (ull i = 0; i + codon_size < domain.size(); i += codon_size)
    {
        string codon = domain.substr(i, codon_size);
        string amino_acid = codon_table.at(codon);
        raw += amino_acid;
    }
    return raw;
}

vector<ui> kmers_encoded(vector<string> kmers)
{
    assert(kmers[0].size() == CasUtil::k);
    vector<ui> encoded;
    for (string kmer : kmers)
    {
        ui kmer_int = 0;
        for (ui i = 0; i < CasUtil::k; i++)
        {
            kmer_int += amino_encoding[kmer[i]] << CasUtil::k * i;
        }
        encoded.push_back(kmer_int);
    }
    return encoded;
}

vector<vector<ull>> cluster_index(const vector<ull>& indices)
{
    vector<vector<ull>> clusters; vector<ull> cluster;
    ull prev = indices[0];
    for (ull index : indices)
    {
        if (index - prev > 5)
        {
            vector<ull> cluster_cp = cluster; clusters.push_back(cluster_cp); cluster.clear();
        }
        cluster.push_back(index); prev = index;
    }
    clusters.push_back(cluster);
    return clusters;
}


const ui cluster_metric_min = 5;

bool good_clusters(const vector<vector<ull>>& clusters)
{
    for (vector<ull> cluster : clusters)
    {   
        if (cluster.size() > cluster_metric_min)
        {
            return true;
        }
    }
    return false;
}

ull demarc_start_clusters(const vector<vector<ull>>& clusters)
{
    for (const vector<ull>& cluster : clusters)
        if (cluster.size() > cluster_metric_min)
            return cluster[0];
    assert(false); return -1;
}

ull demarc_final_clusters(const vector<vector<ull>>& clusters)
{
    for (ull i = clusters.size()-1; i >= 0; i--)
        if (clusters[i].size() > cluster_metric_min)
            return clusters[i][clusters[i].size()-1];
    assert(false); return -1;
}



void compute_demarc(Fragment& frag)
{
    FragDemarc* demarc = new FragDemarc();
    demarc->clust_begin = demarc_start_clusters(frag.clusters);
    demarc->clust_final = demarc_final_clusters(frag.clusters);
    frag.demarc = demarc;
}

void compute_details(Fragment& fragment, const string& genome)
{    
    // ull index_kmer_start = demarc_start_clusters(fragment.clusters);
    // ull index_kmer_end = demarc_final_clusters(fragment.clusters);
    ull index_kmer_start = fragment.demarc->clust_begin;
    ull index_kmer_final = fragment.demarc->clust_final;


    // string protein = fragment.reference_translation->pure.substr(index_kmer_start, (index_kmer_end - index_kmer_start) + CasUtil::k);

    ull raw_pos_start = fragment.reference_translation->pure_mapping[index_kmer_start];
    ull raw_pos_end = fragment.reference_translation->pure_mapping[index_kmer_final];

    ull genome_start;
    ull genome_final;
    if (fragment.reference_translation->pos)
    {
        genome_start = fragment.reference_translation->genome_start + (raw_pos_start * 3); 
        genome_final = fragment.reference_translation->genome_start + ((raw_pos_end + CasUtil::k) * 3) + 3;
    }
    else
    {
        genome_final = fragment.reference_translation->genome_end - (raw_pos_start * 3);
        genome_start = fragment.reference_translation->genome_end - ( ((raw_pos_end + CasUtil::k) * 3) + 3 );
    }
    
    string domain = genome.substr(genome_start, genome_final - genome_start);
    domain = fragment.reference_translation->pos ? domain : Util::reverse_complement(domain);
    string translation = translate_domain(domain);
    string profile = fragment.reference_profile->raw;

    size_t distance = uiLevenshteinDistance(translation, profile);

    FragDetails* details = new FragDetails();
    
    details->genome_start = genome_start;
    details->genome_final = genome_final;
    details->translation = translation;
    details->quality = distance;
    
    fragment.details = details;
}



bool fragment_equivalent(const Fragment& a, const Fragment& b)
{
    return a.demarc->clust_begin == b.demarc->clust_begin && a.demarc->clust_final == b.demarc->clust_final;
}

bool fragment_contains(const Fragment& a, const Fragment& b)
{
    return a.demarc->clust_begin > b.demarc->clust_begin && a.demarc->clust_final < b.demarc->clust_final;
}




vector<Fragment> CasUtil::cas(const vector<CasProfile>& cas_profiles, const vector<Translation>& translations, const string& genome)
{
    fmt::print("\tidentifying cas genes in {} translations...\n", translations.size());
    auto start = time();  

    // ull target_map_size = cas_profiles.size() * translations.size() * 3334; // maximum size of a translation is CasUtil::upstream_size / 3 = 3333.33 (3334), but they are a bit smaller because of stop codons
    // bool* target_map = (bool*) malloc(sizeof(bool) * target_map_size);
    vector<ui> target_map;

    start = time(start, "target map preparation");

    for (ull cas_i = 0; cas_i < cas_profiles.size(); cas_i++)
    {  
        const CasProfile& cas_profile = cas_profiles[cas_i];
        for (ull translation_i = 0; translation_i < translations.size(); translation_i++)
        {
            for (ui i = 0; i < translations[translation_i].pure_kmerized_encoded.size(); i++)
            {
                ui query = translations[translation_i].pure_kmerized_encoded[i];
                auto at_index = cas_profile.hash_table[query % cas_profile.N];
                target_map.push_back(at_index == query);
            }
        }
    }
    start = time(start, "target map construction");

    vector<Fragment> fragments;
    ull target_map_index = 0;
    for (ull cas_i = 0; cas_i < cas_profiles.size(); cas_i++)
    {
        for (ull translation_i = 0; translation_i < translations.size(); translation_i++)
        {
            vector<ull> index;
            for (ull i = 0; i < translations[translation_i].pure_kmerized_encoded.size(); i++)
            {
                if (target_map[target_map_index++])
                {
                    index.push_back(i);
                }
            }

            if (index.size() == 0)
                continue;

            vector<vector<ull>> clusters = cluster_index(index);    
            
            if (!good_clusters(clusters))
                continue;

            Fragment f = {translations[translation_i].reference_crispr, &translations[translation_i], &cas_profiles[cas_i], clusters};
            fragments.push_back(f);
           
        }
    }
    start = time(start, "fragments from target map construction");




    for (Fragment& a : fragments)
        compute_demarc(a);
    start = time(start, "compute fragment demarcs");


    Util::sort(fragments, [](const Fragment& a, const Fragment& b) {
        return a.demarc->clust_begin < b.demarc->clust_begin;
    });
    start = time(start, "cas sort fragments");

    // vector<Fragment> containment_filtered;
    // for (ui i = 0; i < fragments.size(); i++)
    // {
    //     auto a = fragments[i];
    //     bool include_a = true;
    //     for (ui j = 0; j < fragments.size(); j++)
    //     {
    //         if (i == j)
    //             continue;

    //         auto b = fragments[j];
    //         if (fragment_contains(a, b))
    //         {
    //             include_a = false;
    //             break;
    //         }
    //     }

    //     if (include_a)
    //         containment_filtered.push_back(a);
    // }
    // start = time(start, "cas fragment prune");

    for (Fragment& a : fragments)
        compute_details(a, genome);
    start = time(start, "compute fragment details");



    // vector<Fragment> equivalency_filtered;
    // for (ui i = 0; i < containment_filtered.size(); i++)
    // {
    //     auto a = containment_filtered[i];
    //     bool include_a = true;
    //     for (ui j = 0; j < containment_filtered.size(); j++)
    //     {
    //         auto b = containment_filtered[j];
    //         if (i == j) 
    //             continue;
    //         if (fragment_equivalent(a, b))
    //         {
    //             if (a.details->quality < b.details->quality)
    //             {
    //                 include_a = false; break;
    //             }
    //         }
    //     }

    //     if (include_a)
    //         equivalency_filtered.push_back(a);

    // }
    // start = time(start, "cas equivalency prune");


    // Are there multiple fragments aligned to the same name?


    return fragments;
}


void print_fragment(const Fragment& fragment, const string& genome)
{
    fmt::print("\t\t{}\n", fragment.reference_profile->name);
    fmt::print("\t\t\t{}...{}\n", fragment.details->genome_start, fragment.details->genome_final);
    fmt::print("\t\t\t{}\n", fragment.details->quality);
    fmt::print("\t\t\tt:{}\n", fragment.details->translation);
    fmt::print("\t\t\tr:{}\n", fragment.reference_profile->raw);
}

void CasUtil::print_fragments(const vector<Crispr>& crisprs, const vector<Fragment>& fragments, const string& genome)
{

    map<string, vector<Fragment>> crispr_buckets;
    for (const Crispr& c : crisprs)
    {
        string crispr_string = fmt::format("{}:{}\n", c.start, c.k);

        for (const Fragment& f : fragments)
        {
            if (f.reference_crispr->start == c.start && f.reference_crispr->k == c.k)
            {
                crispr_buckets[crispr_string].push_back(f);
            }
        }
    }

    for (auto pairing : crispr_buckets)
    {
        fmt::print("\t{}\n", pairing.first);

        map<string, vector<Fragment>> frag_map;
        for (const Fragment& a : pairing.second)
        {
            frag_map[a.reference_profile->name].push_back(a);
        }

        for (const pair<string, vector<Fragment>>& pairing : frag_map)
        {
            fmt::print("UNIQUE: {}\n", pairing.first);
            for (const Fragment& a : pairing.second)
            {
                print_fragment(a, genome);
            }
        }
    }

}

void CasUtil::write_cache(string file, vector<CasProfile> profiles)
{
    ofstream myfile;
    myfile.open(file);
    for (auto p : profiles)
    {
        myfile << ">" + p.name << endl << p.N << endl;
    }
    myfile.close();
}

map<string, string> loaded_cache;

void CasUtil::load_cache(string file)
{
    loaded_cache = Util::parse_fasta(file);
}

ui CasUtil::get_n(CasProfile& profile)
{
    return atoi(loaded_cache.at(profile.name).c_str());
}

ui CasUtil::gen_n(CasProfile& profile)
{   
    ui N = 1;
    for (;;N++)
    {       
        bool* bools = (bool*) malloc(sizeof(bool) * N);
        memset(bools, 0, sizeof(bool) * N);

        bool succeed = true;
        for (ui kmer : profile.encoded_kmer_set)
        {
            if (bools[kmer % N])
            {
                succeed = false;
                break;
            }
            bools[kmer % N] = true;
        }
        free(bools);
        if (succeed)
            break;    
    }
    return N;
}

vector<CasProfile> prelim_load(string uniprot)
{
    auto start = time();

    auto loading = Util::parse_fasta(uniprot);

    vector<string> names;
    vector<string> raws;
    for (pair<string, string> element : loading) 
    {
        names.push_back(element.first);
        raws.push_back(element.second);
    }

    vector<CasProfile> profiles;
    for (ui i = 0; i < names.size(); i++)
    {
        string name = names[i];
        string raw = raws[i];
        vector<string> kmers = Util::kmerize(raw, CasUtil::k);
        vector<ui> encoded_kmers = kmers_encoded(kmers);
        set<ui> encoded_kmer_set;

        for (ui kmer : encoded_kmers)
        {
            encoded_kmer_set.insert(kmer);
        }

        CasProfile cas_profile
        {
            name,
            raw,
            kmers,
            encoded_kmers,
            encoded_kmer_set,
            nullptr,
            0
        };

        profiles.push_back(cas_profile);

    }
    start = time(start, "crispr profile load");
    return profiles;
}

vector<CasProfile> CasUtil::load(string uniprot, function<ui(CasProfile&)> get_n)
{   
    auto profiles = prelim_load(uniprot);

    for (auto &p : profiles)
    {
        p.N = get_n(p);
        ui* hash_table = (ui*) malloc(sizeof(ui) * p.N);
        memset(hash_table, 0, sizeof(ui) * p.N);
        for (ui kmer : p.encoded_kmer_set)
        {
            hash_table[kmer % p.N] = kmer;
        }
        p.hash_table = hash_table;
    }
    return profiles;
}

void debug_clusters(const vector<vector<ull>>& clusters)
{
    for (vector<ull> cluster : clusters) 
        fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
}


vector<Translation> get_triframe(const string& genome, ull genome_start, ull genome_end, bool pos)
{
    string domain = genome.substr(genome_start, genome_end - genome_start);
    domain = pos ? domain : Util::reverse_complement(domain);
    
    vector<Translation> translations;
    for (ull frame = 0; frame < 3; frame++)
	{
        Translation translation;
        translation.pos = pos;
        translation.genome_start = pos ? genome_start + frame : genome_start;
        translation.genome_end = pos ? genome_end : genome_end - frame;

        translation.raw = translate_domain(domain.substr(frame));

		translation.pure = "";

		ull stop_count = 0;
		ull index = 0;
		for (char elem : translation.raw )
		{
			if (elem == CasUtil::stop_c)
			{
				stop_count++;
				continue;
			}
			translation.pure += elem;
			translation.pure_mapping.push_back(index + stop_count);
			index++;
		}
        translation.pure_kmerized = Util::kmerize(translation.pure, CasUtil::k);
        translation.pure_kmerized_encoded = kmers_encoded(translation.pure_kmerized);
        translations.push_back(translation);
    }
    return translations;
}


vector<Translation> CasUtil::get_translations(const string& genome, const vector<Crispr>& crisprs)
{
    auto start = time();
    vector<Translation> translations;
    for (const Crispr& c : crisprs)
    {
        ull genome_start = c.start - CasUtil::upstream_size; genome_start = genome_start < c.start ? genome_start : 0;
        ull genome_end = c.end + CasUtil::upstream_size; genome_end = genome_end > c.end ? genome_end : genome.size()-1;

        vector<Translation> up = get_triframe(genome, genome_start, c.start, true);
        vector<Translation> down = get_triframe(genome, c.end, genome_end, false);

        for (Translation& t : up)
        {
            t.reference_crispr = &c;
            translations.push_back(t);
        }


        for (Translation& t : down)
        {
            t.reference_crispr = &c;
            translations.push_back(t);
        }
 

    }
    time(start, "get translations");
    return translations;
}
