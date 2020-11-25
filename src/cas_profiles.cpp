#include "cas_profiles.h"

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
	return domain_map[name];
}

vector<CasProfile*> CasProfileUtil::load_profiles(std::filesystem::path directory)
{
	vector<CasProfile*> profiles;
	for (const auto& entry : std::filesystem::directory_iterator(directory))
	{
		string file_path = entry.path().string();
		phmap::BinaryInputArchive archive(file_path.c_str());
		CasProfile* profile = new CasProfile;
		profile->hash_table.load(archive);
		profile->identifier = entry.path().stem().string();
		profiles.push_back(profile);
	}
	fmt::print("loaded {} profiles\n", profiles.size());
	return profiles;
}

// serialization
CasProfile* profile_factory(string id, vector<string> sequences, ll k)
{
	CasProfile* profile = new CasProfile;
	profile->identifier = id;
	for (string& seq : sequences)
	{	
		// clean seq
		seq.erase(remove(seq.begin(), seq.end(), '.'), seq.end());
		seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
		transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

		for (ll i = 0; i < seq.size() - k + 1; i++)
		{
			string kmer = seq.substr(i, k);

			if (kmer.find('O') != string::npos || kmer.find('X') != string::npos || kmer.find('B') != string::npos || kmer.find('Z') != string::npos)
				continue;

			ui kmer_enc = Util::encode_amino_kmer(kmer);
			profile->hash_table.insert(kmer_enc);
		}
	}
	return profile;
}

// serialization
void filter_pfam(std::map<string, string>& domain_map, std::filesystem::path pfam_full, std::filesystem::path pfam_filt) 
{
	ifstream input(pfam_full);
	ofstream output(pfam_filt);
	if (!input.good())
		throw runtime_error(strerror(errno));

	auto non_seq = [](string& line) {
		return
			line.starts_with(" ") ||
			line.starts_with("#=GF") ||
			line.starts_with("#=GS") ||
			line.starts_with("#=GC") ||
			line.starts_with("#=GR") ||
			line.starts_with("# STOCKHOLM 1.0") ||
			line.starts_with("//") ||
			line == "";
	};

	vector<string> seq_buffer;
	string ac;
	ui line_count = 0;
	bool engage = false;
	string line;
	while (getline(input, line))
	{
		if (++line_count % 10000 == 0) fmt::print("{}\n", line_count);

		if (line.starts_with("#=GF AC"))
		{
			ac = line;
			ac.erase(0, 10);
			ac.erase(ac.find_first_of('.'));

			if (domain_map.contains(ac))
			{
				engage = true;
			}

			continue;
		}

		if (!engage)
		{
			continue;
		}

		if (line.starts_with("//"))
		{
			output << ac << endl;

			for (string& line : seq_buffer)
			{
				string sequence = line.substr(line.find_last_of(' ') + 1);

				output << sequence << endl;
			}

			output << "//" << endl;

			seq_buffer.clear();
			engage = false;
			continue;
		}

		if (non_seq(line)) continue;
		seq_buffer.push_back(line);
	}

	input.close();
	output.close();
}

// serialization
vector<CasProfile*> profiles_from_processed_pfam(std::filesystem::path path)
{
	ifstream input(path);

	if (!input.good())
		throw runtime_error(strerror(errno));

	vector<CasProfile*> profiles;
	string line;
	while (getline(input, line))
	{
		string id = line;

		vector<string> sequences;
		while (getline(input, line) && line != "//")
			sequences.push_back(line);

		profiles.push_back(profile_factory(id, sequences, CasProfileUtil::k));
		fmt::print("{} pfam profiles built\n", profiles.size());
	}
	return profiles;
}

// serialization
vector<CasProfile*> profiles_from_tigrfam_dir(std::map<string, string> domain_map, std::filesystem::path tigr_dir)
{
	auto stockholm_to_profile = [](std::filesystem::path stockholm) {
		ifstream file(stockholm);
		if (!file.good())
			throw runtime_error(strerror(errno));

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

		return profile_factory(stockholm.filename().string(), sequences, CasProfileUtil::k);
	};

	vector<CasProfile*> profiles;

	for (const auto& entry : std::filesystem::directory_iterator(tigr_dir))
	{
		if (!domain_map.contains(entry.path().filename().string()))
			continue;

		profiles.push_back(stockholm_to_profile(entry.path()));
		fmt::print("{} tigrfam profiles built\n", profiles.size());
	}
	
	return profiles;
}

// serialization
void CasProfileUtil::serialize()
{	
	std::filesystem::path serialization_dir = "T:/prospector-util/profiles/";
	std::filesystem::path pfam_full = "T:/prospector-util/seed/Pfam-A.full";
	std::filesystem::path pfam_filt = "T:/prospector-util/seed/Pfam-A.filt";
	std::filesystem::path tigrfam_dir = "T:/prospector-util/seed/TIGRFAMs_13.0_SEED";
	std::filesystem::path domain_map_path = "T:/prospector-util/cas/domain_map.tsv";

	// step 0: load domain map
	CasProfileUtil::load_domain_map(domain_map_path);

	// step 1: PFAM full needs to be filtered into PFAM filt
	filter_pfam(domain_map, pfam_full, pfam_filt);
		
	// step 2: generate PFAM profiles
	vector<CasProfile*> profiles_pfam = profiles_from_processed_pfam(pfam_filt);

	// step 3: generate TIGRFAM profiles
	vector<CasProfile*> profiles_tigrfam = profiles_from_tigrfam_dir(domain_map, tigrfam_dir);

	// step 4: serialize all profiles
	auto serialize = [&](CasProfile* profile) {
		std::filesystem::path file_name = serialization_dir / profile->identifier;
		phmap::BinaryOutputArchive archive(file_name.string().c_str());
		profile->hash_table.dump(archive);
	};

	std::for_each(profiles_pfam.begin(), profiles_pfam.end(), serialize);
	std::for_each(profiles_tigrfam.begin(), profiles_tigrfam.end(), serialize);
}


