#include "stdafx.h"
#include "util.h"
#include "cas_profiles.h"
#include <regex>

//void CasUtil::write_cache(string file, vector<CasProfile> profiles)
//{
//    ofstream myfile;
//    myfile.open(file);
//    for (auto p : profiles)
//    {
//        myfile << ">" + p.name << endl << p.N << endl;
//    }
//    myfile.close();
//}

//map<string, string> loaded_cache;
//
//void CasUtil::load_cache(string file)
//{
//    loaded_cache = Util::parse_fasta(file);
//}

//ui CasUtil::get_n(CasProfile& profile)
//{
//    return atoi(loaded_cache.at(profile.name).c_str());
//}

//ui gen_n(const unordered_set<ui>& encoded_kmer_set)
//{
//    ui N = 1;
//    fmt::print("generating N from set size {}... ", encoded_kmer_set.size());
//    for (;; N++)
//    {
//        bool* bools = (bool*)malloc(sizeof(bool) * N);
//        memset(bools, 0, sizeof(bool) * N);
//
//        bool succeed = true;
//        for (ui kmer : encoded_kmer_set)
//        {
//            if (bools[kmer % N])
//            {
//                succeed = false;
//                break;
//            }
//            bools[kmer % N] = true;
//        }
//        free(bools);
//        if (succeed)
//            break;
//    }
//    fmt::print("{}\n", N);
//    return N;
//}

//ui* gen_perfect_hash_table(ui N, unordered_set<ui> encoded_kmer_set)
//{
//    ui* hash_table = (ui*)malloc(sizeof(ui) * N);
//    memset(hash_table, 0, sizeof(ui) * N);
//    for (ui kmer : encoded_kmer_set)
//        hash_table[kmer % N] = kmer;
//    hash_table = hash_table;
//    return hash_table;
//}

//string gn_from_name(const string& name)
//{
//    regex re("cas[0-9]+|csm[0-9]+|csn[0-9]+", regex_constants::icase);
//
//    map<string, string> gn_resolution
//    {
//        {"csn1", "cas9"},
//    };
//
//    cmatch m;
//    bool result = regex_search(name.c_str(), m, re);
//
//    if (!result)
//    {
//        //fmt::print("gn failure on {}\n", name);
//        return "failure";
//    }
//    else
//    {
//        auto final_val = string(m[0]);
//        std::transform(final_val.begin(), final_val.end(), final_val.begin(), ::tolower);
//
//        if (gn_resolution.contains(final_val))
//        {
//            return gn_resolution.at(final_val);
//        }
//        return final_val;
//    }
//}

//void serialize_profile(CasProfile profile)
//{
//    FILE* write_ptr;
//
//    string ting = profile.gn + ".bin";
//    write_ptr = fopen(ting.c_str(), "wb");
//
//    // write first N, then write the entire hashtable
//    ui* n_buffer = (ui*)malloc(1 * sizeof(ui));
//    n_buffer[0] = profile.N;
//    fwrite(n_buffer, sizeof(ui), 1, write_ptr);
//
//    fwrite(profile.hash_table, sizeof(ui), profile.N, write_ptr);
//}

//vector<CasProfile> CasProfileUtil::cas_profiles_from_uniprot_download(string file_path)
//{
//    map<string, string> parsed_fasta = Util::parse_fasta(file_path);
//
//    map<string, unordered_set<ui>> gn_to_encoded_kmers;
//    for (auto const& [name, seq] : parsed_fasta)
//    {
//        string gn = gn_from_name(name);
//
//        if (gn == "failure")
//            continue;
//
//        vector<string> kmers = Util::kmerize(seq, CasProfileUtil::k);
//        for (ui kmer : Util::encode_amino_kmers(kmers, CasProfileUtil::k))
//            gn_to_encoded_kmers[gn].insert(kmer);
//    }
//
//    fmt::print("have gn's:\n");
//    for (auto const& [gn, encoded_kmer_set] : gn_to_encoded_kmers)
//    {
//        fmt::print("{}\n", gn);
//    }
//
//
//    vector<CasProfile> profiles;
//    for (auto const& [gn, encoded_kmer_set] : gn_to_encoded_kmers)
//    {
//        ui N = gen_n(encoded_kmer_set);
//        ui* hash_table = gen_perfect_hash_table(N, encoded_kmer_set);
//        
//        CasProfile profile;
//        profile.gn = gn;
//        profile.N = N;
//        profile.hash_table = hash_table;
//        //profile.hash_table = encoded_kmer_set;
//
//        profiles.push_back(profile);
//    }
//    return profiles;
//}

//vector<CasProfile> CasUtil::load(string uniprot, function<ui(CasProfile&)> get_n)
//{
//    auto start = time();
//
//    auto profiles = prelim_load(uniprot);
//
//    start = time(start, "prelim load");
//
//    for (auto& p : profiles)
//    {
//        p.N = get_n(p);
//        ui* hash_table = (ui*)malloc(sizeof(ui) * p.N);
//        memset(hash_table, 0, sizeof(ui) * p.N);
//        for (ui kmer : p.encoded_kmer_set)
//        {
//            hash_table[kmer % p.N] = kmer;
//        }
//        p.hash_table = hash_table;
//    }
//
//    start = time(start, "postlim load");
//
//    return profiles;
//}

//vector<CasProfile> prelim_load(string uniprot)
//{
//    auto fasta = Util::parse_fasta(uniprot);
//
//    vector<CasProfile> profiles;
//    for (auto pairing : fasta)
//    {
//        string name = pairing.first;
//        string raw = pairing.second;
//        vector<string> kmers = Util::kmerize(raw, k);
//        vector<ui> encoded_kmers = encode_amino_kmers(kmers);
//        set<ui> encoded_kmer_set;
//
//        for (ui kmer : encoded_kmers)
//        {
//            encoded_kmer_set.insert(kmer);
//        }
//
//        string gn = gn_from_name(name);
//
//        if (gn == "failure")
//            continue;
//
//        CasProfile cas_profile
//        {
//            name,
//            gn,
//            raw,
//            kmers,
//            encoded_kmer_set,
//            nullptr,
//            0
//        };
//        profiles.push_back(cas_profile);
//    }
//    return profiles;
//}

//unordered_set<ui> encode_kmer_set(unordered_set<string> kmer_set, ull k)
//{
//	unordered_set<ui> kmer_set_encoded;
//	for (string& kmer : kmer_set)
//	{
//		if (kmer.find('O') != string::npos || kmer.find('X') != string::npos || kmer.find('B') != string::npos || kmer.find('Z') != string::npos)
//			continue;
//
//		kmer_set_encoded.insert(Util::encode_amino_kmer(kmer, k));
//	}
//	return kmer_set_encoded;
//}
//
//unordered_set<string> kmer_set_from_seqs(vector<string> seq_buffer, ull k)
//{
//	unordered_set<string> kmer_set;
//	for (const string& seq : seq_buffer)
//	{
//		for (ull i = 0; i < seq.size() - k + 1; i++)
//		{
//			auto substr = seq.substr(i, k);
//			kmer_set.insert(substr);
//		}
//	}
//	return kmer_set;
//}

void clean_seq(string& seq)
{
	seq.erase(remove(seq.begin(), seq.end(), '.'), seq.end());
	seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

CasProfile* profile_factory(string id, vector<string> seq_buffer, ull k)
{
	CasProfile* profile = new CasProfile;
	profile->gn = id;
	for (string& seq : seq_buffer)
	{	
		for (ull i = 0; i < seq.size() - k + 1; i++)
		{
			string kmer = seq.substr(i, k);

			if (kmer.find('O') != string::npos || kmer.find('X') != string::npos || kmer.find('B') != string::npos || kmer.find('Z') != string::npos)
				continue;

			ui kmer_enc = Util::encode_amino_kmer(kmer);
			//profile->kmer_set.insert(kmer);
			profile->hash_table.insert(kmer_enc);
		}
	}
	return profile;
}

//vector<string> seqs_from_tigrfam_file(string path)
//{
//    ifstream input(path);
//
//    if (!input.good())
//        throw runtime_error(strerror(errno));
//
//    vector<string> seq_buffer;
//
//    string line;
//    while (getline(input, line))
//    {
//        if (line == "" || line.starts_with('#') || line.starts_with('//'))
//            continue;
//
//        auto result = line.find(' ');
//        auto line_without_id = line.substr(result);
//        auto first_not_of = line_without_id.find_first_not_of(' ');
//        auto seq = line_without_id.substr(first_not_of);
//        seq_buffer.push_back(seq);
//    }
//    return seq_buffer;
//}

//vector<const CasProfile*> CasProfileUtil::profiles_from_tigrfam_dir(string dir)
//{
//    map<string, string> map_id_gn
//    {
//		//{"TIGR00287", "cas1"},
//		//{"TIGR03639", "cas1_a"},
//        //{"TIGR01865", "cas9"},
//        //{"TIGR01573", "cas2"},
//        {"TIGR00372", "cas4"},
//    };
//
//    vector<const CasProfile*> profiles;
//    for (const auto& entry : filesystem::directory_iterator(dir))
//    {
//        string id = entry.path().stem().string();
//        if (!map_id_gn.contains(id)) continue;
//		
//		auto seq_buffer = seqs_from_tigrfam_file(entry.path().string());
//
//		auto profile = profile_factory(map_id_gn.at(id), seq_buffer, CasProfileUtil::k);
//
//        profiles.push_back(profile);
//    }
//    return profiles;
//}


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

void CasProfileUtil::pfam_filter(string in, string out)
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
	string id;
	string ac;
	string line;
	ui line_count = 0;
	bool engage = false;
	while (getline(input, line))
	{
		if (++line_count % 10000 == 0) fmt::print("{}\n", line_count);

		if (line.starts_with("#=GF ID"))
		{
			id = line.substr(10);
			continue;
		}

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
			output << id << endl;

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




vector<const CasProfile*> CasProfileUtil::pfam(string path)
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

		//break;
		fmt::print("{} profiles built\n", profiles.size());
	}
	return profiles;
}




void CasProfileUtil::serialize(string dir, vector<const CasProfile*> profiles)
{	
	for (auto profile : profiles)
	{
		string file_name = dir + "\\" + profile->gn;
		phmap::BinaryOutputArchive archive(file_name.c_str());
		profile->hash_table.dump(archive);

	/*	FILE* file = fopen(file_name.c_str(), "wb");
		auto size = profile->hash_table.size();
		vector<ui> buffer;
		for (ui kmer_enc : profile->hash_table)
		{
			buffer.push_back(kmer_enc);
		}

		fwrite(&buffer[0], sizeof(ui), size, file);
	*/
		//fclose(file);
	
	}
}

vector<const CasProfile*> CasProfileUtil::deserialize(string dir)
{

	
	vector<const CasProfile*> profiles;

	for (const auto& entry : filesystem::directory_iterator(dir))
	{
		string file_path = entry.path().string();
		phmap::BinaryInputArchive archive(file_path.c_str());

		CasProfile* profile = new CasProfile;
		profile->hash_table.load(archive);
		profile->gn = entry.path().stem().string();
		profiles.push_back(profile);

		//FILE* file = fopen(file_path.c_str(), "rb");
		//assert(file);
		//fseek(file, 0, SEEK_END);
		//long fsize = ftell(file);
		//fseek(file, 0, SEEK_SET);
		//ui* buffer = (ui*)malloc(fsize);
		//fread(buffer, 4, fsize / 4, file);
		//fclose(file);
		//CasProfile* profile = new CasProfile;
		////profile->hash_table = unordered_set<ui>{ buffer, buffer + (fsize / 4) };
		//profile->hash_table.reserve(fsize / 4);
		//for (ull i = 0; i < fsize / 4; i++)
		//{
		//	profile->hash_table.insert(buffer[i]);
		//}
		//free(buffer);
		//profile->gn = entry.path().stem().string();
		//profiles.push_back(profile);
		//fmt::print("{} profiles deserialized\n", profiles.size());
	
	}
	return profiles;
}

