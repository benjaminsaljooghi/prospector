#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"

#include "fmt/core.h"
#include "fmt/format.h"
#include <bitset>
#include <cstddef>


#define UPSTREAM_SIZE 10000
#define K_FRAGMENT 5


vector<unsigned int> frames{
    0,
    1,
    2
};


const size_t encoding_size = 5;

map<char, int> amino_encoding {
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

int str_to_int(string kmer)
{
    assert(kmer.size() == encoding_size);
    int my_int = 0;
    for (int i = 0; i < encoding_size; i++)
    {
        my_int += amino_encoding[kmer[i]] << encoding_size * i;
    }
    return my_int;
}


vector<int> kmers_encoded(vector<string> kmers)
{
    vector<int> encoded;
    for (string kmer : kmers)
    {
        encoded.push_back(str_to_int(kmer));
    }
    return encoded;
}



class Translation
{
    public:
        size_t genome_start;
        size_t genome_end;
        map<unsigned int, string> translations_raw;
        map<unsigned int, string> translations_pure;
        map<unsigned int, vector<string>> translations_pure_kmerized;
        map<unsigned int, vector<int>> translations_pure_kmerized_encoded;
        map<unsigned int, vector<size_t>> pure_mapping;

        Translation(const string& genome, size_t genome_start, size_t genome_end, unsigned int k, bool rc);
        const char* to_string();
};


Translation::Translation(const string& genome, size_t genome_start, size_t genome_end, unsigned int k, bool rc)

{
    this->genome_start = genome_start;
    this->genome_end = genome_end;
    string domain = genome.substr(this->genome_start, this->genome_end - this->genome_start);
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

    this->translations_raw[0] = seqs[0];
	this->translations_raw[1] = seqs[1];
	this->translations_raw[2] = seqs[2];

	for (auto const& [key, val] : this->translations_raw)
	{
		this->translations_pure[key] = "";
		
		size_t stop_count = 0;
		size_t index = 0;
		for (char elem : val)
		{
			if (elem == STOP_C)
			{
				stop_count++;
				continue;
			}
			this->translations_pure[key] += elem;
			pure_mapping[key].push_back(index + stop_count);
			index++;
		}


        vector<string> kmers = kmerize(this->translations_pure[key], k);
        this->translations_pure_kmerized[key] = kmers;
        this->translations_pure_kmerized_encoded[key] = kmers_encoded(kmers);
	}
}

const char* Translation::to_string()
{
	// size_t frame_count = 3;
	string result = "";
	for (size_t frame : frames)
		result += fmt::format("{}:\n{}\n\n", frame, this->translations_raw[frame]);

	return result.c_str();
}

class CasProfile
{
    public:
        string name;
        string type;
        vector<string> kmers;
        vector<int> encoded_kmers;
		
		CasProfile(string, unsigned int);

        static map<string, vector<CasProfile>> load_casprofiles(string, unsigned int);
};

CasProfile::CasProfile(string _path, unsigned int _k)
{
	this->name = filesystem::path(_path).stem();
	this->kmers = kmerize(parse_fasta_single(_path), _k);
    this->encoded_kmers = kmers_encoded(this->kmers);
    this->type = this->name.substr(0, this->name.find("_"));
}

map<string, vector<CasProfile>> CasProfile::load_casprofiles(string dir, unsigned int k)
{   
    map<string, vector<CasProfile>> profiles;  
    for (const auto& entry : filesystem::directory_iterator(dir))
    {
        CasProfile cas_profile(entry.path(), k);
        profiles[cas_profile.type].push_back(cas_profile);
    }
    return profiles;
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

void print_clusters(const vector<vector<size_t>>& clusters)
{
    for (vector<size_t> cluster : clusters) fmt::print("\t\t {} - {} ({})\n", cluster[0], cluster[cluster.size()-1], cluster.size());
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

struct gene_fragment
{
    const Crispr* reference_crispr;
    const Translation* reference_translation;
    const CasProfile* reference_profile;
    vector<vector<size_t>> clusters;
    size_t frame;
};

bool gene_fragment_less(const gene_fragment& a, const gene_fragment& b)
{
    return demarc_start_clusters(a.clusters) < demarc_start_clusters(b.clusters);
}

size_t gene_fragment_length(const gene_fragment& a)
{
    return demarc_end_clusters(a.clusters) - demarc_start_clusters(a.clusters); 
}

// a is contained within b
bool fragment_contains(const gene_fragment& a, const gene_fragment& b)
{
    size_t a_start = demarc_start_clusters(a.clusters);
    size_t a_end = demarc_end_clusters(a.clusters);
    size_t b_start = demarc_start_clusters(b.clusters);
    size_t b_end = demarc_end_clusters(b.clusters);
    
    bool equivalent = a_start == b_start && a_end == b_end;

    return (!equivalent) && a_start >= b_start && a_end <= b_end;
}

void print_gene_fragment(gene_fragment gene)
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

vector<gene_fragment> detect(const string& genome, const Translation* translation, const CasProfile* cas_profile, const Crispr* crispr)
{
    vector<gene_fragment> fragments;
    for (auto const& [frame, crispr_profile] : translation->translations_pure_kmerized_encoded)
    {
        optional<vector<vector<size_t>>> clusters = detect_single(crispr_profile, cas_profile->encoded_kmers);
        if (clusters)
        {
            gene_fragment fragment = {
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

void cas(const string& genome, const vector<Crispr>& crisprs, string cas_dir)
{
    double start = omp_get_wtime();  

    map<string, vector<CasProfile>> cas_profiles = CasProfile::load_casprofiles(cas_dir, K_FRAGMENT);

    vector<Translation> downstreams;
    vector<Translation> upstreams;

    for (const Crispr& crispr : crisprs)
    {
        // guard against overflow
        size_t genome_start = crispr.start - UPSTREAM_SIZE;
        if (genome_start > crispr.start)
        {
            genome_start = 0;
        }

        Translation down(genome, genome_start, crispr.start, K_FRAGMENT, false);

        size_t genome_end = crispr.end + UPSTREAM_SIZE;
        if (genome_end < crispr.end)
        {
            genome_end = genome.size()-1;
        }
        Translation up(genome, crispr.end, genome_end, K_FRAGMENT, true);
        downstreams.push_back(down);
        upstreams.push_back(up);
    }

    vector<gene_fragment> fragments;

    for (size_t i = 0; i < crisprs.size(); i++)
    {
        for (auto const& [type, profiles] : cas_profiles)
        {
            for (size_t j = 0; j < profiles.size(); j++)
            {
                vector<gene_fragment> __fragments = detect(genome, &downstreams[i], &profiles[j], &crisprs[i]);
                vector<gene_fragment> _fragments = detect(genome, &upstreams[i], &profiles[j], &crisprs[i]);

                fragments.insert(fragments.end(), __fragments.begin(), __fragments.end());
                fragments.insert(fragments.end(), _fragments.begin(), _fragments.end());
            }
        }
    }


    // organize fragments, sorted by position, then sorted by type, then sorted by fragment
    sort(fragments.begin(), fragments.end(), gene_fragment_less);

    // remove any fragments that are a complete subset of any other fragments
    vector<gene_fragment> fragments_filtered;

    for (gene_fragment a : fragments)
    {
        bool do_not_include_a = false;
        for (gene_fragment b : fragments)
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
    

    // place gene fragments into crispr buckets
    for (const Crispr& crispr : crisprs)
    {
        fmt::print("\tcrispr {} {}\n", crispr.start, crispr.k);
        for (const gene_fragment& g : fragments_filtered )
        {
            if (g.reference_crispr->start == crispr.start && g.reference_crispr->k == crispr.k)
            {
                print_gene_fragment(g);
            }
        }
    }


    done(start, "cas");


}

void finish()
{
    fmt::print("terminate called\n");
    exit(0);
}

void debug(vector<Crispr> crisprs, string genome, ui start, ui end)
{

    vector<Crispr> filtered = filter(crisprs, [&](Crispr c) { return c.start > start-100 && c.end < end+100; } );

    sort(filtered.begin(), filtered.end(), CrisprUtil::heuristic_less);

    int how_many = filtered.size();
    for (size_t i = filtered.size()-how_many; i < filtered.size(); i++)
        filtered[i].print(genome);
    finish();
}

vector<string> load_genomes(string dir)
{
    vector<string> genomes;
    for (const auto& entry : filesystem::directory_iterator(dir))
        genomes.push_back(parse_fasta_single(entry.path()));
    return genomes;
}



map<char, ull> _scheme {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}
};
const ull bits = 2;

ull compute(ull e1, ull e2)
{
    ull _xor = (e1 ^ e2);
    ull evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
    ull oddBits = _xor & 0x5555555555555555ull;
    ull comp = (evenBits >> 1) | oddBits;
    return popcount(comp);
}


ull encode(string kmer)
{
    ull encoding = 0;
    for (int i = 0; i < kmer.size(); i++)
        encoding += _scheme[kmer[i]] << i * bits;
    return encoding;
}

void print_encoding(ull encoding)
{
    bitset<64> a(encoding);
    cout << a << endl;
}

int main()
{
    printf("running main...\n");
    double start = omp_get_wtime();

    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";
    const string genome = load_genomes(genome_dir)[0];


    vector<Crispr> crisprs = Prospector::prospector_main(genome);
       
    CrisprUtil::cache_crispr_information(genome, crisprs);

    // debug(crisprs, genome, 1283502, 1283729);

    vector<Crispr> good_heuristic_crisprs = filter(crisprs, [](const Crispr& c) { return c.overall_heuristic >= 0.7; });

    sort(good_heuristic_crisprs.begin(), good_heuristic_crisprs.end(), CrisprUtil::heuristic_greater);

    vector<Crispr> final = CrisprUtil::get_domain_best(good_heuristic_crisprs);

    sort(final.begin(), final.end(), [](const Crispr& a, const Crispr&b) { return a.start < b.start; });

    CrisprUtil::print(genome, final);

  
    cas(genome, final, cas_dir);


    done(start, "main");

    finish();
    return 0;                                                                                                           
}

