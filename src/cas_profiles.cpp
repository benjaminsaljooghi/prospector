#include "cas_profiles.h"

static map<string, string> domain_table;
static vector<CasProfile*> profiles;

void CasProfileUtil::load_profiles(std::filesystem::path path)
{
	for (const auto& entry : std::filesystem::directory_iterator(path))
	{
		string file_path = entry.path().string();
		phmap::BinaryInputArchive archive(file_path.c_str());
		CasProfile* profile = new CasProfile;
		profile->hash_table.load(archive);
		profile->identifier = entry.path().stem().string();
		profiles.push_back(profile);
	}
	fmt::print("loaded {} profiles\n", profiles.size());
}

vector<CasProfile*>& CasProfileUtil::get_profiles()
{
	return profiles;
}

void CasProfileUtil::load_domain_table(std::filesystem::path path)
{
	std::ifstream file(path);
	if (!file)
	{
		throw std::runtime_error("failed to open domain table");
	}

	string line;
	while (getline(file, line))
	{
		auto split = Util::parse(line, "\t");
		domain_table[split[0]] = split[1];
	}
}

bool CasProfileUtil::domain_table_contains(string domain)
{
	return domain_table.contains(domain);
}

string CasProfileUtil::domain_table_fetch(string domain)
{
	return domain_table.at(domain);
}


void clean_seq(string& seq)
{
	seq.erase(remove(seq.begin(), seq.end(), '.'), seq.end());
	seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

CasProfile* profile_factory(string id, vector<string> sequences, ll k)
{
	CasProfile* profile = new CasProfile;
	profile->identifier = id;
	for (string& seq : sequences)
	{	
		clean_seq(seq);

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

void CasProfileUtil::pfam_filter(std::filesystem::path in, std::filesystem::path out, std::filesystem::path pfam_table)
{
	unordered_set<string> known;
	
	std::ifstream list(pfam_table);
	if (!list)
		throw runtime_error(pfam_table.string());

	string line;
	while (getline(list, line))
	{
		known.insert(Util::parse(line, "\t")[0]);
	}
	
	ifstream input(in);
	ofstream output(out);
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
	while (getline(input, line))
	{
		if (++line_count % 10000 == 0) fmt::print("{}\n", line_count);

		if (line.starts_with("#=GF AC"))
		{
			ac = line;
			ac.erase(0, 10);
			ac.erase(ac.find_first_of('.'));

			if (known.contains(ac))
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
				int first_space = line.find(' ');
				line.erase(0, first_space);
				int first_amino = line.find_first_not_of(' ');
				line.erase(0, first_amino);

				clean_seq(line);

				output << line << endl;
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

vector<const CasProfile*> profiles_from_processed_pfam(std::filesystem::path path)
{
	ifstream input(path);

	if (!input.good())
		throw runtime_error(strerror(errno));

	vector<const CasProfile*> profiles;
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


std::map<string, string> parse_tigr_table(std::filesystem::path path)
{
	ifstream file(path);
	if (!file.good())
		throw runtime_error(strerror(errno));

	std::map<string, string> table;

	string line;
	while (getline(file, line))
	{
		auto split = Util::parse(line, "\t");
		table[split[0]] = split[1];
	}

	return table;
}


vector<const CasProfile*> CasProfileUtil::profiles_from_tigrfam_dir(std::filesystem::path tigr_table, std::filesystem::path tigr_dir)
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

			size_t last_space = line.find_last_of(' ');
			string parsable = line.substr(last_space + 1);
			sequences.push_back(parsable);
		}

		return profile_factory(stockholm.filename().string(), sequences, CasProfileUtil::k);
	};

	std::map<string, string> table = parse_tigr_table(tigr_table);

	vector<const CasProfile*> profiles;

	for (const auto& entry : std::filesystem::directory_iterator(tigr_dir))
	{
		if (!table.contains(entry.path().filename().string()))
			continue;

		profiles.push_back(stockholm_to_profile(entry.path()));
		fmt::print("{} tigrfam profiles built\n", profiles.size());
	}
	
	return profiles;
}





void CasProfileUtil::serialize()
{
	std::filesystem::path serialization_dir = "T:/prospector-util/cas/serial_staging/";

	//std::filesystem::path filt = "T:/prospector-util/Pfam-A.filt";
	//vector<const CasProfile*> profiles = profiles_from_processed_pfam(filt);
	
	std::filesystem::path tigrfam_dir = "T:/prospector-util/cas/TIGRFAMs_13.0_SEED";
	std::filesystem::path tigrfam_table = "T:/prospector-util/cas/tigr_to_domain.tsv";

	vector<const CasProfile*> profiles = profiles_from_tigrfam_dir(tigrfam_table, tigrfam_dir);

	for (auto profile : profiles)
	{
		std::filesystem::path file_name = serialization_dir / profile->identifier;
		phmap::BinaryOutputArchive archive(file_name.string().c_str());
		profile->hash_table.dump(archive);
	}
}


