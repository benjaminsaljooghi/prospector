#include "cas_profiles.h"
#include "path.h"
#include "config.h"

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

vector<CasProfile*> generate_pfams(std::filesystem::path dir)
{
	// ifstream input(pfam_full);
	// if (!input.good())
		// throw runtime_error("input not good!");

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

	// vector<string> seq_buffer;
	// string ac;
	// ui line_count = 0;
	// bool engage = false;
	// string line;
	// while (getline(input, line))
	// {
	// 	if (++line_count % 10000 == 0) fmt::print("{}\n", line_count);

	// 	if (line.starts_with("#=GF AC"))
	// 	{
	// 		ac = line;
	// 		ac.erase(0, 10);
	// 		ac.erase(ac.find_first_of('.'));

	// 		cout << ac << endl;
	// 		if (domain_map.contains(ac))
	// 		{
	// 			engage = true;
	// 		}

	// 		continue;
	// 	}

	// 	if (!engage)
	// 	{
	// 		continue;
	// 	}

	// 	if (line.starts_with("//"))
	// 	{
	// 		ofstream output(dir / ac);

	// 		for (string& line : seq_buffer)
	// 		{
	// 			string sequence = line.substr(line.find_last_of(' ') + 1);

	// 			output << sequence << endl;
	// 		}
	// 		output.close();

	// 		seq_buffer.clear();
	// 		engage = false;
	// 		continue;
	// 	}

	// 	if (is_seq(line)) continue;
	// 	seq_buffer.push_back(line);
	// }

	// input.close();

	vector<CasProfile*> profiles;

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
		profiles.push_back(profile_factory(id, sequences, CasProfileUtil::k));
		fmt::print("{} pfam profiles built\n", profiles.size());
	}


	return profiles;
}

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
vector<CasProfile*> generate_from_fasta(std::filesystem::path fasta_dir)
{
	vector<CasProfile*> profiles;
	vector<string> seqs;

	for (const auto& entry : std::filesystem::directory_iterator(fasta_dir))
	{
		string identifier = entry.path().stem().string();

		if (!domain_map.contains(identifier))
		{
			fmt::print("skipping generation of profile {} because the identifier is not in the domain map\n", identifier);
			continue;
		}
		
		std::map<string, string> fasta_seqs = Util::parse_fasta(entry.path().string(), false);

		seqs.clear();
		for (auto const& entry: fasta_seqs) seqs.push_back(entry.second);

		CasProfile* profile = profile_factory(identifier, seqs, CasProfileUtil::k);
		profiles.push_back(profile);
	}

	return profiles;
}

vector<CasProfile*> generate_cogs(std::filesystem::path cog_dir)
{
	return generate_from_fasta(cog_dir);
}

void CasProfileUtil::serialize(std::filesystem::path path_bin_pro)
{	
	auto serialize_profile = [&](CasProfile* profile) {
		std::filesystem::path file_name = path_bin_pro / profile->identifier;
		phmap::BinaryOutputArchive archive(file_name.string().c_str());
		profile->hash_table.dump(archive);
	};

	auto profiles = generate_from_fasta(Config::path_cas_profiles);
	std::for_each(profiles.begin(), profiles.end(), serialize_profile);
}
