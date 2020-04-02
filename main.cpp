#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"

#include "fmt/core.h"
#include "fmt/format.h"
// #include "fmt/format-inl.h"
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
	size_t frame_count = 3;
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

double index_computation = 0;

template <typename T> vector<size_t> build_index_single(const vector<T>& query_kmers, const vector<T>& target_kmers)
{
    vector<size_t> indices;
    double start = omp_get_wtime();
    for (size_t i = 0; i < target_kmers.size(); i++)
    {
        if (contains(query_kmers, target_kmers[i]))
        {
            indices.push_back(i);
        }
    } 
    double end = omp_get_wtime();
    index_computation += end - start;
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
}

size_t demarc_end_clusters(const vector<vector<size_t>>& clusters)
{
    for (size_t i = clusters.size()-1; i >= 0; i--)
        if (clusters[i].size() > 1)
            return clusters[i][clusters[i].size()-1];
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

    // fmt::print("crispr {} {}\n", gene.reference_crispr->start, gene.reference_crispr->k);
    fmt::print("\t\t{}\n", gene.reference_profile->type);
    // fmt::print("\tframe {}\n", gene.frame);
    // fmt::print("\t{}\n", gene.reference_profile->name);
    // fmt::print("\t{} - {}\n", index_kmer_start, index_kmer_end);
    fmt::print("\t\t\t{} - {}\n", genome_start, genome_end);
    fmt::print("\t\t\t{}...{}\n\n", protein.substr(0, 4), protein.substr(protein.length()-4, 4));
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
    map<string, vector<CasProfile>> cas_profiles = CasProfile::load_casprofiles(cas_dir, K_FRAGMENT);

    vector<Translation> downstreams;
    vector<Translation> upstreams;

    for (const Crispr& crispr : crisprs)
    {
        Translation down(genome, crispr.start - UPSTREAM_SIZE, crispr.start, K_FRAGMENT, false);
        Translation up(genome, crispr.end, crispr.end + UPSTREAM_SIZE, K_FRAGMENT, true);
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




    fmt::print("index computation time {}\n", index_computation);
}

void finish()
{
    exit(0);
}

#if DEBUG == 1
void debug(vector<Crispr> crisprs, string genome)
{
    int how_many = crisprs.size();
    for (size_t i = 0; i < how_many; i++)
    {
        if (crisprs[i].k > 29 && crisprs[i].k < 31)
        {
            crisprs[i].print(genome);
        }
    }
    finish();
}
#endif

vector<string> load_genomes(string dir)
{
    vector<string> genomes;
    for (const auto& entry : filesystem::directory_iterator(dir))
        genomes.push_back(parse_fasta_single(entry.path()));
    return genomes;
}





int main()
{

    printf("running main...\n");
    double start = omp_get_wtime();


    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";
    string genome = load_genomes(genome_dir)[1];


    vector<Crispr> crisprs = Prospector::prospector_main(genome);
    
    CrisprUtil::cache_crispr_information(crisprs, genome);
    
    vector<Crispr> good_heuristic_crisprs;
    for (const Crispr& crispr : crisprs)
    {
        #if DEBUG == 0
        if (crispr.overall_heuristic > 0.5)
        #endif
        {
            good_heuristic_crisprs.push_back(crispr);
        }
    }

    sort(good_heuristic_crisprs.begin(), good_heuristic_crisprs.end(), CrisprUtil::heuristic_greater);

    #if DEBUG == 1
        debug(good_heuristic_crisprs, genome);
    #endif

    vector<Crispr> final = CrisprUtil::get_domain_best(good_heuristic_crisprs);

    // map<string, int> spacer_scores = CrisprUtil::get_spacer_scores(domain_best, target_db_path);
    // vector<Crispr> final = CrisprUtil::spacer_score_filtered(domain_best, spacer_scores);
    // CrisprUtil::print(genome, final, spacer_scores);

    sort(final.begin(), final.end(), [](const Crispr& a, const Crispr&b) {
        return a.start < b.start;
    });

    CrisprUtil::print(genome, final);

    double cas_start = omp_get_wtime();    
    cas(genome, final, cas_dir);
    done(cas_start, "cas_gene_detection");

    done(start, "main");

    finish();
    return 0;                                                                                                           
}

