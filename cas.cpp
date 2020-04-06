#include "cas.h"



vector<ui> frames{
    0,
    1,
    2
};


const size_t encoding_size = 5;

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
    {'S', 19},
    {'G', 20},
};

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




ui str_to_int(string kmer)
{
    assert(kmer.size() == encoding_size);
    ui my_int = 0;
    for (ui i = 0; i < encoding_size; i++)
    {
        my_int += amino_encoding[kmer[i]] << encoding_size * i;
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



Translation get_translation(const string& genome, size_t genome_start, size_t genome_end, ui k, bool rc)
{
    string domain = genome.substr(genome_start, genome_end - genome_start);
    domain = rc ? reverse_complement(domain) : domain; 

    size_t codon_size = 3;

    auto frame_shift = [&](string& dna)
    {
        vector<string> amino_acid_seqs{"", "", ""};
        for (size_t frame = 0; frame < 3; frame++)
		{
			for (size_t i = frame; i + codon_size < dna.size(); i += codon_size)
			{
				string codon = dna.substr(i, codon_size);
				string amino_acid = codon_table.at(codon);
				amino_acid_seqs[frame] += amino_acid;
			}
		}
		return amino_acid_seqs;
    };

    vector<string> seqs = frame_shift(domain);

    map<ui, string> translations_raw;
    map<ui, string> translations_pure;
    map<ui, vector<string>> translations_pure_kmerized;
    map<ui, vector<ui>> translations_pure_kmerized_encoded;
    map<ui, vector<size_t>> pure_mapping;

    translations_raw[0] = seqs[0];
	translations_raw[1] = seqs[1];
	translations_raw[2] = seqs[2];

	for (auto const& [key, val] : translations_raw)
	{
		translations_pure[key] = "";
		
		size_t stop_count = 0;
		size_t index = 0;
		for (char elem : val)
		{
			if (elem == STOP_C)
			{
				stop_count++;
				continue;
			}
			translations_pure[key] += elem;
			pure_mapping[key].push_back(index + stop_count);
			index++;
		}


        vector<string> kmers = kmerize(translations_pure[key], k);
        translations_pure_kmerized[key] = kmers;
        translations_pure_kmerized_encoded[key] = kmers_encoded(kmers);
	}

    return Translation
    {
        genome_start,
        genome_end,
        translations_raw,
        translations_pure,
        translations_pure_kmerized,
        translations_pure_kmerized_encoded,
        pure_mapping
    };
}




template <typename T> vector<size_t> build_index_single(const vector<T>& query_kmers, const vector<T>& target_kmers)
{
    vector<size_t> indices;
    for (size_t i = 0; i < target_kmers.size(); i++)
    {
        if (contains(query_kmers, target_kmers[i]))
        {
            indices.push_back(i);
        }
    } 
    return indices;
}

bool good_index(const vector<size_t>& index)
{
    return index.size() > 0;
}

vector<vector<size_t>> cluster_index(const vector<size_t>& indices)
{
    vector<vector<size_t>> clusters; vector<size_t> cluster;
    size_t prev = indices[0];
    for (size_t index : indices)
    {
        if (index - prev > 5)
        {
            vector<size_t> cluster_cp = cluster; clusters.push_back(cluster_cp); cluster.clear();
        }
        cluster.push_back(index); prev = index;
    }
    clusters.push_back(cluster);

    return clusters;
}

bool good_clusters(const vector<vector<size_t>>& clusters)
{
    size_t cluster_requirement = 3;
    bool good_clusters = false;
    for (vector<size_t> cluster : clusters)
    {   
        if (cluster.size() > cluster_requirement)
        {
            good_clusters = true;
            break;
        }
    }
    return good_clusters;
}


template<typename T> optional<vector<vector<size_t>>> detect_single(const vector<T>& crispr_profile, const vector<T>& cas_profile)
{
    vector<size_t> index = build_index_single(cas_profile, crispr_profile);
    if (!good_index(index)) return {};
    vector<vector<size_t>> clusters = cluster_index(index);    
    if (!good_clusters(clusters)) return {};
    return clusters;
}

size_t demarc_start_clusters(const vector<vector<size_t>>& clusters)
{
    for (const vector<size_t>& cluster : clusters)
        if (cluster.size() > 1)
            return cluster[0];
    assert(false); return -1;
}

size_t demarc_end_clusters(const vector<vector<size_t>>& clusters)
{
    for (size_t i = clusters.size()-1; i >= 0; i--)
        if (clusters[i].size() > 1)
            return clusters[i][clusters[i].size()-1];
    assert(false); return(-1);
    return -1;
}


// a is contained within b
bool fragment_contains(const Fragment& a, const Fragment& b)
{
    size_t a_start = demarc_start_clusters(a.clusters);
    size_t a_end = demarc_end_clusters(a.clusters);
    size_t b_start = demarc_start_clusters(b.clusters);
    size_t b_end = demarc_end_clusters(b.clusters);
    
    bool equivalent = a_start == b_start && a_end == b_end;

    return (!equivalent) && a_start >= b_start && a_end <= b_end;
}

void print_gene_fragment(const Fragment& gene)
{
    size_t index_kmer_start = demarc_start_clusters(gene.clusters);
    size_t index_kmer_end = demarc_end_clusters(gene.clusters);
    
    string protein = gene.reference_translation->translations_pure.at(gene.frame).substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);

    size_t raw_pos_start = gene.reference_translation->pure_mapping.at(gene.frame)[index_kmer_start];
    size_t raw_pos_end = gene.reference_translation->pure_mapping.at(gene.frame)[index_kmer_end];

    size_t genome_start = gene.reference_translation->genome_start + (raw_pos_start * 3) + gene.frame;
    size_t genome_end = gene.reference_translation->genome_start + ((raw_pos_end + K_FRAGMENT) * 3) + gene.frame + 3;

    fmt::print("\t\t{} {} - {} {}...{}\n", 
                    gene.reference_profile->type,
                    genome_start,
                    genome_end,
                    protein.substr(0, 4),
                    protein.substr(protein.length()-4, 4)
            );
}

vector<Fragment> detect(const string& genome, const Translation* translation, const CasProfile* cas_profile, const Crispr* crispr)
{
    vector<Fragment> fragments;
    for (auto const& [frame, crispr_profile] : translation->translations_pure_kmerized_encoded)
    {
        optional<vector<vector<size_t>>> clusters = detect_single(crispr_profile, cas_profile->encoded_kmers);
        if (clusters)
        {
            Fragment fragment = {
                crispr,
                translation,
                cas_profile,
                clusters.value(),
                frame
            };
            fragments.push_back(fragment);
        }
    }
    return fragments;
}


vector<Flanks> CasUtil::get_flanks(const string& genome, const vector<Crispr>& crisprs)
{
    auto get_flank = [&](const Crispr& c)
    {
        size_t genome_start = c.start - UPSTREAM_SIZE;
        genome_start = genome_start < c.start ? genome_start : 0;
        Translation up = get_translation(genome, genome_start, c.start, K_FRAGMENT, false);

        size_t genome_end = c.end + UPSTREAM_SIZE;
        genome_end = genome_end > c.end ? genome_end : genome.size()-1;
        Translation down = get_translation(genome, c.end, genome_end, K_FRAGMENT, true);

        Flanks flanks {
            up,
            down,
        };

        return flanks;
    };

    auto start = time();
    vector<Flanks> flanks;
    for (const Crispr& c : crisprs)
        flanks.push_back(get_flank(c));
    time(start, "get flanks");
    return flanks;
}


vector<Fragment> CasUtil::cas(const string& genome, const vector<Crispr>& crisprs, const vector<CasProfile>& cas_profiles, const vector<Flanks>& flanks)
{
    auto start = time();  

    vector<Fragment> fragments;
    for (ull i = 0; i < crisprs.size(); i++)
    {   
        for (ull j = 0; j < cas_profiles.size(); j++)
        {
            vector<Fragment> __fragments = detect(genome, &flanks[i].down, &cas_profiles[j], &crisprs[i]);
            vector<Fragment> _fragments = detect(genome, &flanks[i].up, &cas_profiles[j], &crisprs[i]);
            fragments.insert(fragments.end(), __fragments.begin(), __fragments.end());
            fragments.insert(fragments.end(), _fragments.begin(), _fragments.end());
        }
    }

    // organize fragments, sorted by position, then sorted by type, then sorted by fragment
    sort(fragments.begin(), fragments.end(), [](const Fragment& a, const Fragment& b) {
        return demarc_start_clusters(a.clusters) < demarc_start_clusters(b.clusters);
    });

    // remove any fragments that are a complete subset of any other fragments
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
    
    time(start, "cas");
    return fragments_filtered;
}


void CasUtil::print_fragments(vector<Crispr> crisprs, vector<Fragment> fragments)
{
    for (const Crispr& crispr : crisprs)
    {
        fmt::print("\tcrispr {} {}\n", crispr.start, crispr.k);
        for (const Fragment& g : fragments )
        {
            if (g.reference_crispr->start == crispr.start && g.reference_crispr->k == crispr.k)
            {
                print_gene_fragment(g);
            }
        }
    }
}


vector<CasProfile> CasUtil::load(string dir, ui k)
{   
    vector<CasProfile> profiles;  
    for (const auto& entry : fs::directory_iterator(dir))
    {
        // CasProfile cas_profile(entry.path(), k);
        auto path = entry.path();

        string name = path.stem();
        auto type = name.substr(0, name.find("_"));

        vector<string> kmers = kmerize(parse_fasta_single(path), k);
        vector<ui> encoded_kmers = kmers_encoded(kmers);

        CasProfile cas_profile
        {
            name,
            type,
            kmers,
            encoded_kmers,
        };
        profiles.push_back(cas_profile);
    }
    return profiles;
}



void debug_clusters(const vector<vector<size_t>>& clusters)
{
    for (vector<size_t> cluster : clusters) 
        fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
}