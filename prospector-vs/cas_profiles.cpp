#include "cas_profiles.h"

static map<string, string> domain_table;
static vector<CasProfile*> profiles;

void CasProfileUtil::load_profiles(string path)
{
	for (const auto& entry : filesystem::directory_iterator(path))
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

void CasProfileUtil::load_domain_table(string path)
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

CasProfile* profile_factory(string id, vector<string> seq_buffer, ll k)
{
	CasProfile* profile = new CasProfile;
	profile->identifier = id;
	for (string& seq : seq_buffer)
	{	
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

bool non_seq(string& line)
{
	return
		line.starts_with(" ") ||
		line.starts_with("#=GF") ||
		line.starts_with("#=GS") ||
		line.starts_with("#=GC") ||
		line.starts_with("#=GR") ||
		line.starts_with("# STOCKHOLM 1.0") ||
		line.starts_with("//") ||
		line == "";
}

void pfam_filter(string in, string out)
{
	unordered_set<string> known
	{
		"PF18765",
		"PF18759",
		"PF18541",
		"PF18525",
		"PF18516",
		"PF18510",
		"PF18501",
		"PF18470",
		"PF18395",
		"PF18320",
		"PF18211",
		"PF18186",
		"PF18182",
		"PF18070",
		"PF18061",
		"PF18019",
		"PF17955",
		"PF17953",
		"PF17952",
		"PF17894",
		"PF17893",
		"PF17783",
		"PF17667",
		"PF17289",
		"PF17262",
		"PF17167",
		"PF16813",
		"PF16595",
		"PF16593",
		"PF16592",
		"PF14301",
		"PF14279",
		"PF14239",
		"PF13936",
		"PF13715",
		"PF13707",
		"PF13671",
		"PF13532",
		"PF13424",
		"PF13416",
		"PF13412",
		"PF13395",
		"PF13384",
		"PF13366",
		"PF13307",
		"PF13280",
		"PF13245",
		"PF13180",
		"PF13175",
		"PF13087",
		"PF13086",
		"PF12900",
		"PF12848",
		"PF12840",
		"PF12802",
		"PF12705",
		"PF12469",
		"PF11327",
		"PF10902",
		"PF10040",
		"PF09827",
		"PF09822",
		"PF09711",
		"PF09709",
		"PF09707",
		"PF09706",
		"PF09704",
		"PF09703",
		"PF09702",
		"PF09701",
		"PF09700",
		"PF09670",
		"PF09659",
		"PF09658",
		"PF09657",
		"PF09652",
		"PF09651",
		"PF09623",
		"PF09620",
		"PF09618",
		"PF09617",
		"PF09615",
		"PF09614",
		"PF09611",
		"PF09609",
		"PF09559",
		"PF09485",
		"PF09484",
		"PF09481",
		"PF09455",
		"PF09369",
		"PF09344",
		"PF09339",
		"PF09002",
		"PF08843",
		"PF08798",
		"PF08696",
		"PF08448",
		"PF08388",
		"PF08281",
		"PF08279",
		"PF08223",
		"PF08220",
		"PF07848",
		"PF07717",
		"PF07715",
		"PF07118",
		"PF06884",
		"PF06733",
		"PF06209",
		"PF06165",
		"PF06023",
		"PF05635",
		"PF05598",
		"PF05107",
		"PF05065",
		"PF04851",
		"PF04828",
		"PF04464",
		"PF04446",
		"PF04434",
		"PF04408",
		"PF04198",
		"PF04140",
		"PF03787",
		"PF03750",
		"PF03703",
		"PF03601",
		"PF03237",
		"PF02954",
		"PF02796",
		"PF02661",
		"PF02562",
		"PF02518",
		"PF02342",
		"PF02321",
		"PF01978",
		"PF01966",
		"PF01930",
		"PF01905",
		"PF01881",
		"PF01867",
		"PF01844",
		"PF01797",
		"PF01710",
		"PF01638",
		"PF01636",
		"PF01627",
		"PF01614",
		"PF01609",
		"PF01584",
		"PF01547",
		"PF01047",
		"PF01032",
		"PF01022",
		"PF00929",
		"PF00575",
		"PF00535",
		"PF00455",
		"PF00313",
		"PF00271",
		"PF00270",
		"PF00221",
		"PF00196",
		"PF00158",
		"PF00117",
		"PF00080",
		"PF00078",
		"PF00072",
		"PF00005",
	};
	
	ifstream input(in);
	ofstream output(out);
	if (!input.good())
		throw runtime_error(strerror(errno));

	vector<string> seq_buffer;
	//string id;
	string ac;
	string line;
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


vector<const CasProfile*> profiles_from_processed_pfam(string path)
{
	ifstream input(path);

	if (!input.good())
		throw runtime_error(strerror(errno));

	vector<const CasProfile*> profiles;
	string line;
	while (getline(input, line))
	{
		string id = line;

		vector<string> seq_buffer;
		while (getline(input, line) && line != "//")
		{
			seq_buffer.push_back(line);
		}

		profiles.push_back(profile_factory(id, seq_buffer, CasProfileUtil::k));
		fmt::print("{} profiles built\n", profiles.size());
	}
	return profiles;
}




void CasProfileUtil::serialize()
{
	string in = "T:\\crispr\\cas\\Pfam-A.full";
	string out = "T:\\crispr\\cas\\Pfam-A.filt";
	pfam_filter(in, out);
	vector<const CasProfile*> profiles = profiles_from_processed_pfam(out);
	
	for (auto profile : profiles)
	{
		string file_name = "T:\\crispr\\cas\\serial_staging\\" + profile->identifier;
		phmap::BinaryOutputArchive archive(file_name.c_str());
		profile->hash_table.dump(archive);
	}
	
	exit(0);
}