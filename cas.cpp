#include "cas.h"

#define STOP "Z"
#define STOP_C 'Z'
#define ENCODING_SIZE 5

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
    {"TAA", STOP},
    {"TAG", STOP},
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
    {"TGA", STOP},
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

ui str_to_int(string kmer)
{
    assert(kmer.size() == ENCODING_SIZE);
    ui my_int = 0;
    for (ui i = 0; i < ENCODING_SIZE; i++)
    {
        my_int += amino_encoding[kmer[i]] << ENCODING_SIZE * i;
    }
    return my_int;
}

vector<ui> kmers_encoded(vector<string> kmers)
{
    vector<ui> encoded;
    for (string kmer : kmers)
    {
        encoded.push_back(str_to_int(kmer));
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

bool good_clusters(const vector<vector<ull>>& clusters)
{
    ull cluster_requirement = 3;
    bool good_clusters = false;
    for (vector<ull> cluster : clusters)
    {   
        if (cluster.size() > cluster_requirement)
        {
            good_clusters = true;
            break;
        }
    }
    return good_clusters;
}

ull demarc_start_clusters(const vector<vector<ull>>& clusters)
{
    for (const vector<ull>& cluster : clusters)
        if (cluster.size() > 1)
            return cluster[0];
    assert(false); return -1;
}

ull demarc_end_clusters(const vector<vector<ull>>& clusters)
{
    for (ull i = clusters.size()-1; i >= 0; i--)
        if (clusters[i].size() > 1)
            return clusters[i][clusters[i].size()-1];
    assert(false); return(-1);
    return -1;
}


// a is contained within b
bool fragment_contains(const Fragment& a, const Fragment& b)
{
    ull a_start = demarc_start_clusters(a.clusters);
    ull a_end = demarc_end_clusters(a.clusters);
    ull b_start = demarc_start_clusters(b.clusters);
    ull b_end = demarc_end_clusters(b.clusters);
    
    bool equivalent = a_start == b_start && a_end == b_end;

    return (!equivalent) && a_start >= b_start && a_end <= b_end;
}


// struct CasSet
// {
//     vector<ui> profiles;
//     vector<ui> starts;
//     vector<ui> sizes;
//     ui profile_count;
// };

// struct CrisprSet
// {
//     vector<ui> profiles;
//     vector<ui> starts;
//     vector<ui> sizes;
//     ui profile_count;

// };

// CrisprSet crisprset_from_translations(const vector<Translation>& translations)
// {
//     CrisprSet crisprset;
//     ui _crispr_profile_start = 0;
//     for (ull translation_i = 0; translation_i < translations.size(); translation_i++)
//     {
//         auto profile = translations[translation_i].pure_kmerized_encoded;
//         crisprset.starts.push_back(_crispr_profile_start);
//         crisprset.profiles.insert(crisprset.profiles.end(), profile.begin(), profile.end());
//         crisprset.sizes.push_back(profile.size());
//         _crispr_profile_start += profile.size();
//     }
//     return crisprset;
// }

// CasSet casset_from_profiles(const vector<CasProfile>& cas_profiles)
// {
//     CasSet casset;
//     ui _cas_profile_start = 0;
//     casset.profile_count = cas_profiles.size();
//     for (ull cas_i = 0; cas_i < casset.profile_count; cas_i++)
//     {
//         auto cas_profile = cas_profiles[cas_i].encoded_kmers;
//         casset.starts.push_back(_cas_profile_start);
//         casset.profiles.insert(casset.profiles.end(), cas_profile.begin(), cas_profile.end());
//         casset.sizes.push_back(cas_profile.size());
//         _cas_profile_start += cas_profile.size();
//     }
//     return casset;
// }

vector<Fragment> CasUtil::cas(const vector<CasProfile>& cas_profiles, const vector<Translation>& translations)
{
    fmt::print("\tidentifying cas genes in {} translations...\n", translations.size());
    auto start = time();  

    // ull target_map_size = cas_profiles.size() * translations.size() * 3334; // maximum size of a translation is UPSTREAM_SIZE / 3 = 3333.33 (3334), but they are a bit smaller because of stop codons
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

    Util::sort(fragments, [](const Fragment& a, const Fragment& b) {
        return demarc_start_clusters(a.clusters) < demarc_start_clusters(b.clusters);
    });

    start = time(start, "cas sort fragments");
    vector<Fragment> fragments_filtered;
    for (Fragment a : fragments)
    {
        bool do_not_include_a = false;
        for (Fragment b : fragments)
        {
            if (fragment_contains(a, b))
            {
                do_not_include_a = true; break;
            }
        }

        if (!do_not_include_a)
        {
            fragments_filtered.push_back(a);
        }
    }
    start = time(start, "cas fragment prune");

    return fragments_filtered;
}


void print_fragment(const Fragment& fragment, const string& genome)
{
    ull index_kmer_start = demarc_start_clusters(fragment.clusters);
    ull index_kmer_end = demarc_end_clusters(fragment.clusters);
    
    // string protein = fragment.reference_triframe->translations[fragment.frame].pure.substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);

    // ull raw_pos_start = fragment.reference_triframe->translations[fragment.frame].pure_mapping[index_kmer_start];
    // ull raw_pos_end = fragment.reference_triframe->translations[fragment.frame].pure_mapping[index_kmer_end];

    // ull genome_start = fragment.reference_triframe->genome_start + (raw_pos_start * 3) + fragment.frame;
    // ull genome_end = fragment.reference_triframe->genome_start + ((raw_pos_end + K_FRAGMENT) * 3) + fragment.frame + 3;

    string protein = fragment.reference_translation->pure.substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);

    ull raw_pos_start = fragment.reference_translation->pure_mapping[index_kmer_start];
    ull raw_pos_end = fragment.reference_translation->pure_mapping[index_kmer_end];

    ull genome_start;
    ull genome_end;
    if (fragment.reference_translation->pos)
    {
        genome_start = fragment.reference_translation->genome_start + (raw_pos_start * 3); 
        genome_end = fragment.reference_translation->genome_start + ((raw_pos_end + K_FRAGMENT) * 3) + 3;
    }
    else
    {
        genome_end = fragment.reference_translation->genome_end - (raw_pos_start * 3);
        genome_start = fragment.reference_translation->genome_end - ( ((raw_pos_end + K_FRAGMENT) * 3) + 3 );
    }
    

    string pro_begin = protein.substr(0, 4);
    string pro_end = protein.substr(protein.length()-4);

    string domain = genome.substr(genome_start, genome_end - genome_start);
    domain = fragment.reference_translation->pos ? domain : Util::reverse_complement(domain);
    string trans = translate_domain(domain);
    string trans_begin = trans.substr(0, 4);
    string trans_end = trans.substr(trans.length()-4);

    string reference_begin = fragment.reference_profile->raw.substr(0, 4);
    string reference_end= fragment.reference_profile->raw.substr(fragment.reference_profile->raw.length()-4);

    fmt::print("\t\t{} {} {} {} - {} {}...{} {}...{} {}...{}\n", 
                    (reference_begin == trans_begin && reference_end == trans_end) ? "[ OK ]" : "[ BAD ]",
                    (pro_begin == trans_begin && pro_end == trans_end) ? "[ OK ]" : "[ BAD ]",
                    fragment.reference_profile->name,
                    genome_start,
                    genome_end,
                    pro_begin,
                    pro_end,
                    trans_begin,
                    trans_end,
                    reference_begin,
                    reference_end
            );
}

void CasUtil::print_fragments(const vector<Crispr>& crisprs, const vector<Fragment>& fragments, const string& genome)
{
    for (const Crispr& crispr : crisprs)
    {
        fmt::print("\tcrispr {} {}\n", crispr.start, crispr.k);
        for (const Fragment& g : fragments )
        {
            if (g.reference_crispr->start == crispr.start && g.reference_crispr->k == crispr.k)
            {
                print_fragment(g, genome);
            }
        }
    }
}


vector<CasProfile> CasUtil::load(string dir, ui k)
{   
    auto start = time();

    vector<CasProfile> profiles;  
    for (const auto& entry : fs::directory_iterator(dir))
    {
        // CasProfile cas_profile(entry.path(), k);
        auto path = entry.path();

        string name = path.stem();
        auto type = name.substr(0, name.find("_"));

        string raw = Util::parse_fasta_single(path);
        vector<string> kmers = Util::kmerize(raw, k);
        vector<ui> encoded_kmers = kmers_encoded(kmers);
        set<ui> encoded_kmer_set;

        for (ui kmer : encoded_kmers)
        {
            encoded_kmer_set.insert(kmer);
        }

        ui N = 1;
        for (;;N++)
        {       
            bool* bools = (bool*) malloc(sizeof(bool) * N);
            memset(bools, 0, sizeof(bool) * N);
   
            bool succeed = true;
            for (ui kmer : encoded_kmer_set)
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

        ui* hash_table = (ui*) malloc(sizeof(ui) * N);
        memset(hash_table, 0, sizeof(ui) * N);
        for (ui kmer : encoded_kmer_set)
        {
            hash_table[kmer % N] = kmer;
        }

        CasProfile cas_profile
        {
            name,
            type,
            raw,
            kmers,
            encoded_kmers,
            encoded_kmer_set,
            hash_table,
            N
        };
        profiles.push_back(cas_profile);
    }
    start = time(start, "crispr profile load");
    return profiles;
}



void debug_clusters(const vector<vector<ull>>& clusters)
{
    for (vector<ull> cluster : clusters) 
        fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
}


vector<Translation> get_triframe(const string& genome, ull genome_start, ull genome_end, ui k, bool pos)
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
			if (elem == STOP_C)
			{
				stop_count++;
				continue;
			}
			translation.pure += elem;
			translation.pure_mapping.push_back(index + stop_count);
			index++;
		}
        translation.pure_kmerized = Util::kmerize(translation.pure, k);
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
        ull genome_start = c.start - UPSTREAM_SIZE; genome_start = genome_start < c.start ? genome_start : 0;
        ull genome_end = c.end + UPSTREAM_SIZE; genome_end = genome_end > c.end ? genome_end : genome.size()-1;

        vector<Translation> up = get_triframe(genome, genome_start, c.start, K_FRAGMENT, true);
        vector<Translation> down = get_triframe(genome, c.end, genome_end, K_FRAGMENT, false);

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
