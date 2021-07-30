#include "cas_profiles.h"
#include "path.h"

static std::map<string, string> domain_map;

void CasProfileUtil::load_domain_map(std::filesystem::path path)
{
	std::ifstream file(path);
	if (!file)
		throw std::runtime_error("failed to open domain table");

	string line;
	while (getline(file, line))
	{
		auto split = Util::parse(line, "\t");
		domain_map[split[0]] = split[1];
	}
}

bool CasProfileUtil::domain_table_contains(string name)
{
	return domain_map.contains(name);
}

string CasProfileUtil::domain_table_fetch(string name)
{
	std::map<string, string>::iterator it = domain_map.find(name);
	if (it == domain_map.end())
	{
		fmt::print("UNKNOWN IDENTIFIER: {}\n", name);
		throw new exception();
	}
	return it->second;
}

// bool CasProfileUtil::domain_contained(string query_domain)
// {
// 	for (auto const& [identifer, domain] : domain_map)
// 		if (domain.find(query_domain) != string::npos)
// 			return true;
// 	return false;
// }

std::map<string, string> CasProfileUtil::get_domain_map()
{
	return domain_map;
}

void CasProfileUtil::print_profiles(vector<CasProfile*> profiles)
{
    for (CasProfile* p : profiles)
    {
        string id = p->identifier;
        string domain = CasProfileUtil::domain_table_fetch(p->identifier);
        fmt::print("{}\t{}\n", id, domain);
    }   
}

vector<CasProfile*> CasProfileUtil::deserialize_profiles(std::filesystem::path directory)
{
	vector<CasProfile*> profiles;
	for (const auto& entry : std::filesystem::directory_iterator(directory))
	{
		string file_path = entry.path().string();
		phmap::BinaryInputArchive archive(file_path.c_str());
		CasProfile* profile = new CasProfile;
		profile->hash_table.load(archive);
		profile->identifier = entry.path().stem().string();

		if (!CasProfileUtil::domain_table_contains(profile->identifier))
		{
			fmt::print("deserialization of {} skipped due to domain table containment failure\n", profile->identifier);
			continue;
		}
			
		profiles.push_back(profile);
	}
	CasProfileUtil::print_profiles(profiles);

	fmt::print("loaded {} profiles\n", profiles.size());
	return profiles;
}

CasProfile* profile_factory(string id, vector<string> sequences, ull k)
{
	CasProfile* profile = new CasProfile;
	profile->identifier = id;
	for (string& seq : sequences)
	{	
		// clean seq
		seq.erase(remove(seq.begin(), seq.end(), '.'), seq.end());
		seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
		transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

		for (ull i = 0; i < seq.size() - k + 1; i++)
		{
			string kmer = seq.substr(i, k);

			if (kmer.find('O') != string::npos || kmer.find('X') != string::npos || kmer.find('B') != string::npos || kmer.find('Z') != string::npos)
				continue;

			try
			{
				ui kmer_enc = Util::encode_amino_kmer(kmer);
				profile->hash_table.insert(kmer_enc);
			}
			catch(const std::exception& e)
			{
				std::cerr << e.what() << '\n';
			}
		}
	}

	double sum = 0;
	for (ui i = 0; i < sequences.size(); i++)
		sum += sequences[i].size();

	
	std::sort(sequences.begin(), sequences.end(), [] (const std::string& first, const std::string& second){ return first.size() < second.size(); });
	profile->length_median = sequences[sequences.size() / 2].length();
	profile->length_mean = sum / sequences.size();
	profile->length_min = sequences[0].size();
	profile->length_max = sequences[sequences.size()-1].size();

	fmt::print("{}\t{}\t{}\t{}\t{}\n", id, profile->length_median, profile->length_mean, profile->length_min, profile->length_max);

	return profile;
}

vector<CasProfile*> get_profiles_pfam(map<string, vector<string>>& pfam_id_to_seqs) 
{
	vector<CasProfile*> profiles;
	for (auto const& [pfam_id, seqs] : pfam_id_to_seqs)
	{
		profiles.push_back(profile_factory(pfam_id, seqs, CasProfileUtil::k));
	}
	return profiles;
}

map<string, vector<string>> get_seqs_pfam(std::filesystem::path dir)
{
	auto is_seq = [](string& line) {
		return !(
			line.starts_with(" ") ||
			line.starts_with("#=GF") ||
			line.starts_with("#=GS") ||
			line.starts_with("#=GC") ||
			line.starts_with("#=GR") ||
			line.starts_with("# STOCKHOLM 1.0") ||
			line.starts_with("//") ||
			line == "");
	};

	map<string, vector<string>> id_to_seqs; // this is like PF1234: ">a CGACAGTCAG     >b ACAGCATCGACATCGACA" etc (except they are amino seqs)
    for (const auto& entry : std::filesystem::directory_iterator(dir))
	{
		ifstream raw(entry.path().string());

		if (!raw.good())
			throw runtime_error("input not good!");

		string line;
		string id = entry.path().stem().string();
		vector<string> sequences;

		while (getline(raw, line))
		{
			if (is_seq(line))
			{
				string sequence = line.substr(line.find_last_of(' ') + 1);
				sequences.push_back(sequence);
			}
		}
		id_to_seqs[id] = sequences; 
	}
	return id_to_seqs;
}



// map<string, vector<string>> get_seqs_from_fasta_dir(std::filesystem::path dir)
// {
// 	map<string, vector<string>> alles_seqs; // remember, we are parsing a fasta, but the specific ids in each fasta file do not matter, we want to link the FILE'S NAME (COG123) with a vector of all sequences in the file.
// 	for (const auto& entry : std::filesystem::directory_iterator(dir))
// 	{
// 		string identifier = entry.path().stem().string();
// 		std::map<string, string> fasta_seqs = Util::parse_fasta(entry.path().string(), false);
// 		for (auto const& fast_seq: fasta_seqs)
// 			alles_seqs[identifier].push_back(fast_seq.second);
// 	}
// 	return alles_seqs;
// }

// vector<string> get_seqs_tigrfam(std::filesystem::path dir)
// {
// 	daigdad;
// }

vector<CasProfile*> generate_tigrfams(std::filesystem::path tigr_dir)
{
	auto stockholm_to_profile = [](std::filesystem::path stockholm, string identifier) {
		ifstream file(stockholm);
		if (!file.good())
			throw runtime_error("input not good!");

		std::unordered_set<char> skip{ '#', '/' };

		vector<string> sequences;
		string line;
		while (getline(file, line))
		{
			if (skip.contains(line[0]))
				continue;

			string sequence = line.substr(line.find_last_of(' ') + 1);
			sequences.push_back(sequence);
		}

		return profile_factory(identifier, sequences, CasProfileUtil::k);
	};

	vector<CasProfile*> profiles;

	for (const auto& entry : std::filesystem::directory_iterator(tigr_dir))
	{

		string identifier = entry.path().filename().string();

		if (identifier.ends_with(".SEED"))
			identifier = identifier.substr(0, 9);

		if (!domain_map.contains(identifier))
			continue;

		CasProfile* profile = stockholm_to_profile(entry.path(), identifier);
		profiles.push_back(profile);
		fmt::print("built: {}\n", profile->identifier);
		fmt::print("{} tigrfam profiles built\n", profiles.size());
	}
	return profiles;
}

// This needs to be versatile in that it can handle fasta seqs that are both single-line and multi-line
// vector<CasProfile*> generate_from_fasta(std::filesystem::path fasta_dir)
// {
	
// }

// vector<CasProfile*> generate_cogs(std::filesystem::path cog_dir)
// {
// 	return generate_from_fasta(cog_dir);
// }

// std::filesystem::path path_bin_pro = "/home/ben/crispr/prospector-data/bin_pro/"; // duplicated in main
std::filesystem::path makarova_dir = "/home/ben/crispr/prospector-data/raw_pro/makarova/";
// std::filesystem::path pfam_full = data_root / "seed/Pfam-A.full/Pfam-A.full";
std::filesystem::path pfam_dir = "/home/ben/crispr/prospector-data/raw_pro/PFAMS/";
// std::filesystem::path tigrfam_dir = data_root / "seed/TIGRFAMs_14.0_SEED";
std::filesystem::path cog_dir = "/home/ben/crispr/prospector-data/raw_pro/COG/";

// vector<CasProfile*> generate_profiles_from_master_list_of_all_seqs()
// {
// 	if (!domain_map.contains(identifier))
// 	{
// 		fmt::print("skipping generation of profile {} because the identifier is not in the domain map\n", identifier);
// 		continue;
// 	}
// }

void CasProfileUtil::serialize(std::filesystem::path path_bin_pro)
{	
	auto serialize_profile = [&](CasProfile* profile) {
		std::filesystem::path file_name = path_bin_pro / profile->identifier;
		phmap::BinaryOutputArchive archive(file_name.string().c_str());
		profile->hash_table.dump(archive);
	};


	// get all seqs
	map<string, vector<string>> pfam_seqs = get_seqs_pfam(pfam_dir);
	// map<string, vector<string>> cog_seqs = get_seqs_from_fasta_dir(cog_dir);
	// map<string, vector<string>> makarova_seqs = get_seqs_from_fasta_dir(makarova_dir);
	// map<string, vector<string>> tigrfram_seqs = get_seqs_tigrfam(tigrfam_dir);

	// generate_profiles_from_master_list_of_all_seqs();



	// auto profiles_makarova = generate_from_fasta(makarova_dir);
	// auto profiles_pfam = generate_pfams(pfam_dir);
	auto profiles_pfam = get_profiles_pfam(pfam_seqs);
	// auto profiles_tigrfams = generate_tigrfams(tigrfam_dir);
	// auto profiles_cog = generate_cogs(cog_dir);

	// std::for_each(profiles_makarova.begin(), profiles_makarova.end(), serialize_profile);
	std::for_each(profiles_pfam.begin(), profiles_pfam.end(), serialize_profile);
	// std::for_each(profiles_tigrfams.begin(), profiles_tigrfams.end(), serialize_profile);
	// std::for_each(profiles_cog.begin(), profiles_cog.end(), serialize_profile);
}
