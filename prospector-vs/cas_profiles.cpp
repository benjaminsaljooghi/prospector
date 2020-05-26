#include "stdafx.h"
//#include "crispr.h"
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

ui gen_n(const unordered_set<ui>& encoded_kmer_set)
{
    ui N = 1;
    fmt::print("generating N from set size {}... ", encoded_kmer_set.size());
    for (;; N++)
    {
        bool* bools = (bool*)malloc(sizeof(bool) * N);
        memset(bools, 0, sizeof(bool) * N);

        bool succeed = true;
        for (ui kmer : encoded_kmer_set)
        {
            if (bools[kmer % N])
            {
                succeed = false;
                break;
            }
            bools[kmer % N] = true;
        }
        free(bools);
        if (succeed)
            break;
    }
    fmt::print("{}\n", N);
    return N;
}

ui* gen_perfect_hash_table(ui N, unordered_set<ui> encoded_kmer_set)
{
    ui* hash_table = (ui*)malloc(sizeof(ui) * N);
    memset(hash_table, 0, sizeof(ui) * N);
    for (ui kmer : encoded_kmer_set)
        hash_table[kmer % N] = kmer;
    hash_table = hash_table;
    return hash_table;
}

string gn_from_name(const string& name)
{
    regex re("cas[0-9]+|csm[0-9]+|csn[0-9]+", regex_constants::icase);

    map<string, string> gn_resolution
    {
        {"csn1", "cas9"},
    };

    cmatch m;
    bool result = regex_search(name.c_str(), m, re);

    if (!result)
    {
        //fmt::print("gn failure on {}\n", name);
        return "failure";
    }
    else
    {
        auto final_val = string(m[0]);
        std::transform(final_val.begin(), final_val.end(), final_val.begin(), ::tolower);

        if (gn_resolution.contains(final_val))
        {
            return gn_resolution.at(final_val);
        }
        return final_val;
    }
}

vector<string> multi_seq_from_tigfram(string file_path)
{
    ifstream input(file_path);

    if (!input.good())
    {
        throw runtime_error(strerror(errno));
    }

    vector<string> seqs;

    string line;
    while (getline(input, line))
    {
        if (line == "" || line.starts_with('#') || line.starts_with('//'))
        {
            continue;
        }

        auto result = line.find(' ');
        auto line_without_id = line.substr(result);
        auto first_not_of = line_without_id.find_first_not_of(' ');
        auto seq = line_without_id.substr(first_not_of);
        seq.erase(std::remove(seq.begin(), seq.end(), '.'), seq.end());
        seqs.push_back(seq);
    }

    return seqs;
}


unordered_set<string> kmer_set_from_tigfram(string file_path, ull k)
{
    auto seqs = multi_seq_from_tigfram(file_path);
    unordered_set<string> kmer_set;
    for (const string& seq : seqs)
    {
        for (ull i = 0; i < seq.size() - k + 1; i++)
        {
            auto substr = seq.substr(i, k);
            kmer_set.insert(substr);
        }
    }
    return kmer_set;
}

vector<CasProfile> CasProfileUtil::cas_profiles_from_tigrfam(string dir)
{
    map<string, string> map_id_gn
    {
        //{"TIGR01865", "cas9"},
        {"TIGR01573", "cas2"},
        //{"TIGR00372", "cas4"},
    };

    vector<CasProfile> profiles;

    for (const auto& entry : filesystem::directory_iterator(dir))
    {
        string id = entry.path().stem().string();

        if (!map_id_gn.contains(id)) continue;

        unordered_set<string> kmer_set = kmer_set_from_tigfram(entry.path().string(), CasProfileUtil::k);

        unordered_set<ui> hash_table;
        for (string kmer : kmer_set)
            hash_table.insert(Util::encode_amino_kmer(kmer, CasProfileUtil::k));

        CasProfile profile;
        profile.gn = map_id_gn.at(id);
        profile.hash_table = hash_table;
        profile.kmer_set = kmer_set;

        //profile.N = gen_n(hash_table);
        //profile.hash_table = gen_perfect_hash_table(profile.N, hash_table);
        //profile.hash_table = hash_table;
        //profile.hash_table = kmer_set;

        profiles.push_back(profile);
    }
        
    return profiles;
}

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




