#include "cas.h"



vector<ui> frames{
    0,
    1,
    2
};


const ull encoding_size = 5;

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



TriFrame get_triframe(const string& genome, ull genome_start, ull genome_end, ui k, bool rc)
{
    TriFrame triframe;
    triframe.genome_start = genome_start;
    triframe.genome_end = genome_end;

    string domain = genome.substr(genome_start, genome_end - genome_start);
    domain = rc ? reverse_complement(domain) : domain; 

    ull codon_size = 3;

    auto frame_shift = [&](string& dna)
    {
        vector<string> amino_acid_seqs{"", "", ""};
        for (ull frame = 0; frame < 3; frame++)
		{
			for (ull i = frame; i + codon_size < dna.size(); i += codon_size)
			{
				string codon = dna.substr(i, codon_size);
				string amino_acid = codon_table.at(codon);
				amino_acid_seqs[frame] += amino_acid;
			}
		}
		return amino_acid_seqs;
    };

    vector<string> raws = frame_shift(domain);

	
    for (string raw : raws)
	{
        Translation translation;

        translation.raw = raw;
		translation.pure = "";

		ull stop_count = 0;
		ull index = 0;
		for (char elem : raw)
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
        translation.pure_kmerized = kmerize(translation.pure, k);
        translation.pure_kmerized_encoded = kmers_encoded(translation.pure_kmerized);
        triframe.translations.push_back(translation);
    }
    return triframe;
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

void print_fragment(const Fragment& fragment)
{
    ull index_kmer_start = demarc_start_clusters(fragment.clusters);
    ull index_kmer_end = demarc_end_clusters(fragment.clusters);
    
    string protein = fragment.reference_triframe->translations[fragment.frame].pure.substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);

    ull raw_pos_start = fragment.reference_triframe->translations[fragment.frame].pure_mapping[index_kmer_start];
    ull raw_pos_end = fragment.reference_triframe->translations[fragment.frame].pure_mapping[index_kmer_end];

    ull genome_start = fragment.reference_triframe->genome_start + (raw_pos_start * 3) + fragment.frame;
    ull genome_end = fragment.reference_triframe->genome_start + ((raw_pos_end + K_FRAGMENT) * 3) + fragment.frame + 3;

    fmt::print("\t\t{} {} - {} {}...{}\n", 
                    fragment.reference_profile->name,
                    genome_start,
                    genome_end,
                    protein.substr(0, 4),
                    protein.substr(protein.length()-4, 4)
            );
}

vector<Fragment> CasUtil::cas(const string& genome, const vector<Crispr>& crisprs, const vector<CasProfile>& cas_profiles, const vector<Flanks>& flanks)
{
    fmt::print("\tidentifying cas genes in {} crisprs...\n", crisprs.size());
    auto start = time();  



    // ull num_runs = cas_profiles.size() * crisprs.size() * 3; // running each cas profile against each of the three triframes per crispr. Will need to be expanded to six frames later.

    vector<ui> _crispr_profiles;
    vector<ui> _crispr_profile_starts;
    vector<ui> _crispr_profile_sizes;
    ui _crispr_profile_start = 0;
    ui _crispr_profile_count = crisprs.size() * 3;
    ui _crispr_count = crisprs.size();
    for (ull crispr_i = 0; crispr_i < _crispr_count; crispr_i++)
    {
        for (ull frame = 0; frame < 3; frame++)
        {
            auto crispr_profile = flanks[crispr_i].up.translations[frame].pure_kmerized_encoded;

            _crispr_profile_starts.push_back(_crispr_profile_start);
            _crispr_profiles.insert(_crispr_profiles.end(), crispr_profile.begin(), crispr_profile.end());
            _crispr_profile_sizes.push_back(crispr_profile.size());
            _crispr_profile_start += crispr_profile.size();
        }
    }

    vector<ui> _cas_profiles;
    vector<ui> _cas_profile_starts;
    vector<ui> _cas_profile_sizes;
    ui _cas_profile_start = 0;
    ui _cas_profile_count = cas_profiles.size();
    for (ull cas_i = 0; cas_i < _cas_profile_count; cas_i++)
    {
        auto cas_profile = cas_profiles[cas_i].encoded_kmers;

        _cas_profile_starts.push_back(_cas_profile_start);
        _cas_profiles.insert(_cas_profiles.end(), cas_profile.begin(), cas_profile.end());
        _cas_profile_sizes.push_back(cas_profile.size());
        _cas_profile_start += cas_profile.size();
    }

    ull target_map_size = cas_profiles.size() * crisprs.size() * 3 * 3334; // 3 = frames, and maximum size of a crispr is UPSTREAM_SIZE / 3 = 3333.33 (3334)
    bool* target_map = (bool*) malloc(sizeof(bool) * target_map_size);

    start = time(start, "target map preparation");
    ull target_map_index = 0;
    for (ull cas_i = 0; cas_i < _cas_profile_count; cas_i++)
    {
        for (ull crispr_i = 0; crispr_i < _crispr_count; crispr_i++)
        {
            for (ull frame = 0; frame < 3; frame++)
            {
                ui start = _crispr_profile_starts[crispr_i * 3 + frame];
                ui size = _crispr_profile_sizes[crispr_i * 3 + frame];
                for (ui i = 0; i < size; i++)
                {
                    target_map[target_map_index++] = contains(&_cas_profiles[0] + _cas_profile_starts[cas_i], _cas_profile_sizes[cas_i], _crispr_profiles[start + i]);
                }
            }
        }
    }
    start = time(start, "target map construction");

    
    target_map_index = 0;
    vector<Fragment> all_fragments;
    for (ull cas_i = 0; cas_i < _cas_profile_count; cas_i++)
    {
        for (ull crispr_i = 0; crispr_i < _crispr_count; crispr_i++)
        {
            for (ull frame = 0; frame < 3; frame++)
            {
                ui size = _crispr_profile_sizes[crispr_i * 3 + frame];
                vector<ull> index;
                for (ull i = 0; i < size; i++)
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

                Fragment fragment = {
                    &(crisprs[crispr_i]),
                    &(flanks[crispr_i].up),
                    &(cas_profiles[cas_i]),
                    clusters,
                    frame
                };
                all_fragments.push_back(fragment);
            }
           
        }
    }

    // target_map_index = 0;
    // ull frame = 0;
    // for (ull cas_i = 0; cas_i < cas_profiles.size(); cas_i++)
    // {
    //     for (ull crispr_i = 0; crispr_i < crisprs.size(); crispr_i++)
    //     {
    //         for (ull frame = 0; frame < 3; frame++)
    //         {
    //             // auto crispr_profile = flanks[crispr_i].up.translations[frame].pure_kmerized_encoded;
    //             vector<ull> index;    
    //             for (ull i = 0; i < _crispr_profile_sizes[]; i++)
    //             {
    //                 if (target_map[target_map_index++])
    //                 {
    //                     index.push_back(i);
    //                 }
    //             }

    //             if (index.size() == 0)
    //                 continue;

    //             vector<vector<ull>> clusters = cluster_index(index);    
                
    //             if (!good_clusters(clusters))
    //                 continue;

    //             Fragment fragment = {
    //                 &(crisprs[crispr_i]),
    //                 &(flanks[crispr_i].up),
    //                 &(cas_profiles[cas_i]),
    //                 clusters,
    //                 0
    //                 // frame++ % 3
    //             };
    //             all_fragments.push_back(fragment);

    //         }
    //     }
    // }

    start = time(start, "fragments from target map construction");

    sort(all_fragments, [](const Fragment& a, const Fragment& b) {
        return demarc_start_clusters(a.clusters) < demarc_start_clusters(b.clusters);
    });

    start = time(start, "cas sort fragments");
    vector<Fragment> fragments_filtered;
    for (Fragment a : all_fragments)
    {
        bool do_not_include_a = false;
        for (Fragment b : all_fragments)
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


void CasUtil::print_fragments(vector<Crispr> crisprs, vector<Fragment> fragments)
{
    for (const Crispr& crispr : crisprs)
    {
        fmt::print("\tcrispr {} {}\n", crispr.start, crispr.k);
        for (const Fragment& g : fragments )
        {
            if (g.reference_crispr->start == crispr.start && g.reference_crispr->k == crispr.k)
            {
                print_fragment(g);
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



void debug_clusters(const vector<vector<ull>>& clusters)
{
    for (vector<ull> cluster : clusters) 
        fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
}


vector<Flanks> CasUtil::get_flanks(const string& genome, const vector<Crispr>& crisprs)
{
    auto get_flank = [&](const Crispr& c)
    {
        ull genome_start = c.start - UPSTREAM_SIZE;
        genome_start = genome_start < c.start ? genome_start : 0;

        ull genome_end = c.end + UPSTREAM_SIZE;
        genome_end = genome_end > c.end ? genome_end : genome.size()-1;

        TriFrame up = get_triframe(genome, genome_start, c.start, K_FRAGMENT, false);
        TriFrame down = get_triframe(genome, c.end, genome_end, K_FRAGMENT, true);

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
