#include "stdafx.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "cas.h"
#include "time.h"
#include "debug.h"

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

ui gen_n(const set<ui>& encoded_kmer_set)
{
    ui N = 1;
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
    return N;
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
        fmt::print("gn failure on {}\n", name);
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
        if (line.starts_with('#') || line.starts_with('//'))
        {
            continue;
        }


        // get first space
        auto result = line.find(' ');

        auto line_without_id = line.substr(result);

        auto first_not_of = line_without_id.find_first_not_of(' ');


        auto seq = line_without_id.substr(first_not_of);
        // then continue until next char is not a space

        seqs.push_back(seq);
    }

    return seqs;
}

vector<string> kmerize_tigfram_seq(string seq, ull k)
{
    vector<string> kmers;
    for (ull i = 0; i < seq.size() - k + 1; i++)
    {
        auto substr = seq.substr(i, k);

        if (substr.find('.') != string::npos)
        {
            continue;
        }

        kmers.push_back(substr);
    }
    return kmers;
}

vector<string> kmer_set_from_tigfram(string file_path, ull k)
{
    auto seqs = multi_seq_from_tigfram(file_path);
    vector<string> kmer_set;
    for (const string& seq : seqs)
    {
        auto kmers = kmerize_tigfram_seq(seq, k);
        kmer_set.insert(kmer_set.end(), kmers.begin(), kmers.end());
    }
    return kmer_set;
}

void serialize_profile(CasProfile profile)
{
    FILE* write_ptr;

    string ting = profile.gn + ".bin";
    write_ptr = fopen(ting.c_str(), "wb");

    // write first N, then write the entire hashtable
    ui* n_buffer = (ui*)malloc(1 * sizeof(ui));
    n_buffer[0] = profile.N;
    fwrite(n_buffer, sizeof(ui), 1, write_ptr);


    fwrite(profile.hash_table, sizeof(ui), profile.N, write_ptr);
}


vector<ui> encode_amino_kmers(vector<string> kmers)
{
    assert(kmers[0].size() == CasUtil::k);
    vector<ui> encoded(kmers.size());
    memset(&encoded[0], 0, sizeof(ui) * kmers.size());
    for (ui j = 0; j < kmers.size(); j++)
    {
        string kmer = kmers[j];
        for (ui i = 0; i < CasUtil::k; i++)
        {
            encoded[j] += Util::amino_encoding.at(kmer[i]) << CasUtil::k * i;
        }
    }
    return encoded;
}

void cas_seqs_from_uniprot_download(string file_path)
{
    map<string, string> parsed_fasta = Util::parse_fasta(file_path);

    map<string, set<ui>> gn_to_encoded_kmers;
    for (auto const& [name, seq] : parsed_fasta)
    {
        string gn = gn_from_name(name);
        vector<string> kmers = Util::kmerize(seq, CasUtil::k);
        for (ui kmer : encode_amino_kmers(kmers))
        {
            gn_to_encoded_kmers[gn].insert(kmer);
        }
    }

    // now generate the profiles
    vector<CasProfile> profiles;

    for (auto const& [gn, kmer_set] : gn_to_encoded_kmers)
    {
        ui N = gen_n(kmer_set);

        // generate profile
    }


}




vector<CasProfile> prelim_load(string uniprot)
{
    auto fasta = Util::parse_fasta(uniprot);

    vector<CasProfile> profiles;
    for (auto pairing : fasta)
    {
        string name = pairing.first;
        string raw = pairing.second;
        vector<string> kmers = Util::kmerize(raw, CasUtil::k);
        vector<ui> encoded_kmers = kmers_encoded(kmers);
        set<ui> encoded_kmer_set;

        for (ui kmer : encoded_kmers)
        {
            encoded_kmer_set.insert(kmer);
        }

        string gn = gn_from_name(name);

        if (gn == "failure")
            continue;

        CasProfile cas_profile
        {
            name,
            gn,
            raw,
            kmers,
            encoded_kmer_set,
            nullptr,
            0
        };
        profiles.push_back(cas_profile);
    }
    return profiles;
}

vector<CasProfile> CasUtil::load(string uniprot, function<ui(CasProfile&)> get_n)
{
    auto start = time();

    auto profiles = prelim_load(uniprot);

    start = time(start, "prelim load");

    for (auto& p : profiles)
    {
        p.N = get_n(p);
        ui* hash_table = (ui*)malloc(sizeof(ui) * p.N);
        memset(hash_table, 0, sizeof(ui) * p.N);
        for (ui kmer : p.encoded_kmer_set)
        {
            hash_table[kmer % p.N] = kmer;
        }
        p.hash_table = hash_table;
    }

    start = time(start, "postlim load");

    return profiles;
}



