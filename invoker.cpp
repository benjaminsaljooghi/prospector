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



class CrisprProfile
{
    public:

        const Crispr& crispr;
        const Translation translation;
        CrisprProfile(Crispr& _crispr, Translation _translation);
};


CrisprProfile::CrisprProfile(Crispr& _crispr, Translation _translation)
:	crispr(_crispr),
	translation(_translation)
{
}





class ProfileExecution
{
    public:
        CasProfile& cas_profile;
        CrisprProfile& crispr_profile;
        vector<size_t> index;
        vector<vector<size_t>> index_clusters;
        bool good_execution = false;
        unsigned int good_frame = 0;

        ProfileExecution(CasProfile& _cas_profile, CrisprProfile& _crispr_profile)
        :	cas_profile(_cas_profile),
            crispr_profile(_crispr_profile)
        {}

        void build_index()
        {
            for (unsigned int frame : frames)
            {
                vector<string> crispr_kmers = this->crispr_profile.translation.translations_pure_kmerized.at(frame);
                vector<string> cas_kmers = this->cas_profile.kmers;

                vector<size_t> indices;
                printf("comparing %zd kmers against %zd kmers\n", crispr_kmers.size(), cas_kmers.size());
                for (size_t i = 0; i < crispr_kmers.size(); i++)
                {
                    if ( contains(cas_kmers, crispr_kmers[i]) )
                    {
                        indices.push_back(i);
                    }
                } 

                if (indices.size() == 0)
                    continue;

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

                if (good_clusters)
                {
                    this->good_execution = true;
                    this->good_frame = frame;
                    this->index = indices;
                    this->index_clusters = clusters;
                    printf("good index found at frame %d\n", frame);
                    string a = this->crispr_profile.translation.to_string().c_str();
                    printf("%s\n", );
                    break;
                }

            }

        }

        void interpret(string genome)
        {
            assert(this->good_execution);

            printf("profile %s; CRISPR %d %d\n", cas_profile.name.c_str(), crispr_profile.crispr.start, crispr_profile.crispr.k);

            printf("\tframe: %d\n", this->good_frame);

            vector<vector<size_t>> clusters = this->index_clusters;

            // underlying cluster information
            // for (vector<size_t> cluster : clusters)
                // printf("\t \t %zd - %zd (%zd)\n", cluster[0], cluster[cluster.size()-1], cluster.size());

            vector<size_t>  demarc_start;
            vector<size_t>  demarc_end;
            for (vector<size_t> cluster : clusters)
            {
                if (cluster.size() > 1)
                {
                    demarc_start = cluster; break;
                }
            }

            for (size_t i = clusters.size()-1; i >= 0; i--)
            {
                if (clusters[i].size() > 1)
                {
                    demarc_end = clusters[i]; break;
                }
            }

            size_t index_kmer_start = demarc_start[0];
            size_t index_kmer_end = demarc_end[demarc_end.size()-1];
            string demarcated_amino = this->crispr_profile.translation.translations_pure.at(this->good_frame).substr(index_kmer_start, (index_kmer_end - index_kmer_start) + K_FRAGMENT);



            size_t raw_pos_start = this->crispr_profile.translation.pure_mapping.at(this->good_frame)[index_kmer_start];
            size_t raw_pos_end = this->crispr_profile.translation.pure_mapping.at(this->good_frame)[index_kmer_end];


            size_t genome_upstream_start = this->crispr_profile.crispr.start - UPSTREAM_SIZE; // WRONG. WILL NOT WORK FOR RC.
            size_t genome_start = genome_upstream_start + (raw_pos_start * 3) + this->good_frame;
            size_t genome_end = genome_upstream_start + ((raw_pos_end + K_FRAGMENT) * 3) + this->good_frame + 3; // not sure why I need this final +3


            printf("\t \t %zd -  %zd (%zd - %zd) \n", index_kmer_start, index_kmer_end, genome_start, genome_end );
            printf("%s\n", demarcated_amino.c_str());

        }

};






void cas(vector<CrisprProfile> crispr_profiles, map<string, vector<CasProfile>> cas_profile_set, string genome)
{
    double start = omp_get_wtime();

    map<string, vector<ProfileExecution>> executions;

    vector<ProfileExecution> pending_interpretation;

    for (CrisprProfile crispr_profile : crispr_profiles)
    {
        for (auto const& [type, cas_profiles] : cas_profile_set)
        {
            for (CasProfile cas_profile : cas_profiles)
            {
                ProfileExecution execution = ProfileExecution(cas_profile, crispr_profile);
                execution.build_index();
                if (execution.good_execution)
                {
                    // the CasProfile of type TYPE is good enough for this Crispr, therefore no need to check
                    // any more CasProfiles of this type for THIS Crispr.
                    pending_interpretation.push_back(execution);
                    break;
                }
            }
        }
    }


    for (ProfileExecution execution : pending_interpretation)
    {
        execution.interpret(genome);
    }


    done(start, "cas detection");
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


    vector<CrisprProfile> crispr_profiles;
    for (Crispr& crispr : final)
    {
        string region = genome.substr(crispr.start - UPSTREAM_SIZE, UPSTREAM_SIZE);
        Translation translation(region, K_FRAGMENT);
        CrisprProfile crispr_profile(crispr, translation);
        crispr_profiles.push_back(crispr_profile);
    }


    cas(crispr_profiles, cas_profiles, genome);
    
    // try inverse
    crispr_profiles.clear();

    for (Crispr& crispr : final)
    {
        string region = genome.substr(crispr.end, UPSTREAM_SIZE);
        region = reverse_complement(region);
        Translation translation(region, K_FRAGMENT);
        CrisprProfile crispr_profile(crispr, translation);
        crispr_profiles.push_back(crispr_profile);
    }

    cas(crispr_profiles, cas_profiles, genome);
    


    done(start, "invoker");

    finish();
    return 0;                                                                                                           
}

