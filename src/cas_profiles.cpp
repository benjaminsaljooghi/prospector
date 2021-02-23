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
	return domain_map[name];
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
		profiles.push_back(profile);
	}
	for (size_t i = 0; i < profiles.size(); i++)
		fmt::print("loaded: {}\n", profiles[i]->identifier);

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

			ui kmer_enc = Util::encode_amino_kmer(kmer);
			profile->hash_table.insert(kmer_enc);
		}
	}
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

vector<CasProfile*> generate_cogs(std::filesystem::path cog_dir)
{
	vector<CasProfile*> profiles;

	for (const auto& entry : std::filesystem::directory_iterator(cog_dir))
	{
		string identifier = entry.path().stem().string();

		if (!domain_map.contains(identifier))
			continue;

		std::map<string, string> fasta_seqs = Util::parse_fasta(entry.path().string());
		vector<string> seqs;

		for (auto const& entry : fasta_seqs)
			seqs.push_back(entry.second);

		CasProfile* profile = profile_factory(identifier, seqs, CasProfileUtil::k);
		profiles.push_back(profile);
	}

	return profiles;
}


std::filesystem::path data_root = "/home/ben/crispr/data/";

std::filesystem::path serialization_dir = data_root / "profiles/"; // duplicated in main
std::filesystem::path pfam_full = data_root / "seed/Pfam-A.full/Pfam-A.full";
std::filesystem::path pfam_dir = data_root / "seed/PFAMS";
std::filesystem::path tigrfam_dir = data_root / "seed/TIGRFAMs_14.0_SEED";
std::filesystem::path cog_dir = data_root / "seed/cog/";
void CasProfileUtil::serialize()
{	
	auto serialize_profile = [&](CasProfile* profile) {
		std::filesystem::path file_name = serialization_dir / profile->identifier;
		phmap::BinaryOutputArchive archive(file_name.string().c_str());
		profile->hash_table.dump(archive);
	};

	auto profiles_pfam = generate_pfams(pfam_dir);
	// auto profiles_tigrfams = generate_tigrfams(tigrfam_dir);
	// auto profiles_cog = generate_cogs(cog_dir);

	std::for_each(profiles_pfam.begin(), profiles_pfam.end(), serialize_profile);
	// std::for_each(profiles_tigrfams.begin(), profiles_tigrfams.end(), serialize_profile);
	// std::for_each(profiles_cog.begin(), profiles_cog.end(), serialize_profile);
}
