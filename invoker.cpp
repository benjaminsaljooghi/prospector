#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "blast.h"

#include "fmt/core.h"
#include "fmt/format.h"
// #include "fmt/format-inl.h"

#define UPSTREAM_SIZE 10000
#define K_FRAGMENT 5


vector<unsigned int> frames{
    0,
    1,
    2
};




class Translation
{
    public:
        string& nucleotide_sequence;
        map<unsigned int, string> translations_raw;
        // map<string, vector<string>> translations_raw_kmerized;
        map<unsigned int, string> translations_pure;
        map<unsigned int, vector<string>> translations_pure_kmerized;
        map<unsigned int, vector<size_t>> pure_mapping;

        Translation(string&, unsigned int k);
         
        string to_string();
};


Translation::Translation(string& _seq, unsigned int k)
:
	nucleotide_sequence(_seq)
{
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

    vector<string> seqs = frame_shift(nucleotide_sequence);

    this->translations_raw[0] = seqs[0];
	this->translations_raw[1] = seqs[1];
	this->translations_raw[2] = seqs[2];

	// the genome positions are cached and computed here, they are not computed on the fly
	// they are cached in the object via a kind of "map"
	// that is, if I ask for the "index" of a translation,
	// it gets first mapped to the index in the translation_raw
	// and then that index is multiplied by 3 to get me 
	// the genome index

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


		this->translations_pure_kmerized[key] = kmerize(this->translations_pure[key], k);
	}


}


string Translation::to_string()
{
	size_t frame_count = 3;

	string result = "";

	for (size_t frame = 0; frame < frame_count; frame++)
	{
		result += fmt::format("{}:\n{}\n\n", frame, this->translations_raw[frame]);
	}

	return result;
}



class CasProfile
{
    public:
        string name;
        string type;
        vector<string> kmers;
		
		CasProfile(string, unsigned int);

        static map<string, vector<CasProfile>> load_casprofiles(string, unsigned int);
};

CasProfile::CasProfile(string _path, unsigned int _k)
{
	this->name = filesystem::path(_path).stem();
	this->kmers = kmerize(parse_fasta_single(_path), _k);
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


vector<size_t> build_index_single(const vector<string>& query_kmers, const vector<string>& target_kmers)
{
    vector<size_t> indices;
    // printf("comparing %zd kmers against %zd kmers\n", query_kmers.size(), target_kmers.size());
    for (size_t i = 0; i < target_kmers.size(); i++)
    {
        if ( contains(query_kmers, target_kmers[i]) )
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
    // underlying cluster information
    for (vector<size_t> cluster : clusters)
        printf("\t \t %zd - %zd (%zd)\n", cluster[0], cluster[cluster.size()-1], cluster.size());
}

optional<vector<vector<size_t>>> detect_single(const vector<string>& crispr_profile, const vector<string>& cas_profile)
{
    vector<size_t> index = build_index_single(cas_profile, crispr_profile);

    if (!good_index(index))
        return {};

    vector<vector<size_t>> clusters = cluster_index(index);    

    if (!good_clusters(clusters))
        return {};

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

bool detect(const string& genome, const Translation& translation, string cas_profile_name, const vector<string>& cas_profile, const Crispr& crispr)
{
    for (auto const& [frame, crispr_profile] : translation.translations_pure_kmerized)
    {
        optional<vector<vector<size_t>>> clusters = detect_single(crispr_profile, cas_profile);
        if (clusters)
        {
            size_t index_kmer_start = demarc_start_clusters(clusters.value());
            size_t index_kmer_end = demarc_end_clusters(clusters.value());
            
            printf("profile %s; CRISPR %d %d\n", cas_profile_name.c_str(), crispr.start, crispr.k);
            printf("\tframe: %d\n", frame);

            string demarcated_amino = translation.translations_pure.at(frame).substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);

            size_t raw_pos_start = translation.pure_mapping.at(frame)[index_kmer_start];
            size_t raw_pos_end = translation.pure_mapping.at(frame)[index_kmer_end];

            size_t genome_upstream_start = crispr.start - UPSTREAM_SIZE; // WRONG. WILL NOT WORK FOR RC.
            size_t genome_start = genome_upstream_start + (raw_pos_start * 3) + frame;
            size_t genome_end = genome_upstream_start + ((raw_pos_end + K_FRAGMENT) * 3) + frame + 3;

            printf("\t \t %zd -  %zd (%zd - %zd) \n", index_kmer_start, index_kmer_end, genome_start, genome_end );
            printf("%s\n", demarcated_amino.c_str());
            return true;
        }
    }
    return false;
}

void cas(const string& genome, const vector<Crispr>& crisprs, const map<string, vector<CasProfile>>& cas_profiles )
{
    size_t num = crisprs.size();

    vector<Translation> downstreams;
    vector<Translation> upstreams;

    for (size_t i = 0; i < num; i++)
    {
        string down = genome.substr(crisprs[i].start - UPSTREAM_SIZE, UPSTREAM_SIZE);
        string up = genome.substr(crisprs[i].end, UPSTREAM_SIZE);
        downstreams.push_back(Translation(down, K_FRAGMENT));
        upstreams.push_back(Translation(up, K_FRAGMENT));
    }

    for (size_t i = 0; i < num; i++)
    {
        for (auto const& [type, profiles] : cas_profiles)
        {
            for (CasProfile profile : profiles)
            {
                bool good = detect(genome, downstreams[i], profile.name, profile.kmers, crisprs[i]);
                if (good) break;
            }
        }
    }
}

void finish()
{
    exit(0);
}

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




vector<string> load_genomes(string dir)
{

    vector<string> genomes;
    for (const auto& entry : filesystem::directory_iterator(dir))
        genomes.push_back(parse_fasta_single(entry.path()));
    return genomes;
}


int main()
{


    printf("running invoker...\n");
    double start = omp_get_wtime();

    string genome_dir = "crispr-data/genome";
    string cas_dir = "crispr-data/cas";
    string target_db_path = "crispr-data/phage/bacteriophages.fasta";
    string genome = load_genomes(genome_dir)[1];



    // int __start = 1830078;
    // int __end = 1833441;
    // string domain = genome.substr(__start, __end-__start);
    // domain = reverse_complement(domain);
    // Translation domain_translation(domain, K_FRAGMENT);
    // printf("%s\n", domain_translation.to_string().c_str());
    // finish();


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

    CrisprUtil::print(genome, final);

    // cas
    map<string, vector<CasProfile>> cas_profiles = CasProfile::load_casprofiles(cas_dir, K_FRAGMENT);

    cas(genome, final, cas_profiles);

    done(start, "invoker");

    finish();
    return 0;                                                                                                           
}

