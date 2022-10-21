#include "cas_profiles.h"
#include "path.h"
#include "config.h"
#include <boost/algorithm/string.hpp>

static std::map<string, string> domain_map;

void CasProfileUtil::load_domain_map(std::filesystem::path path) {
    std::ifstream file(path);
    if (!file)
        throw std::runtime_error("failed to open domain table");

    string line;
    while (getline(file, line)) {
        auto split = Util::parse(line, "\t");
        domain_map[split[0]] = split[1];
    }
}

bool CasProfileUtil::domain_table_contains(string name) {
    return domain_map.contains(name);
}

string CasProfileUtil::domain_table_fetch(string name) {
    std::map<string, string>::iterator it = domain_map.find(name);
    if (it == domain_map.end()) {
        fmt::print("UNKNOWN IDENTIFIER: {}\n", name);
        throw new exception();
    }
    return it->second;
}

std::map<string, string> CasProfileUtil::get_domain_map() {
    return domain_map;
}

void CasProfileUtil::print_profiles(vector<CasProfile *> profiles) {
    for (CasProfile *p: profiles) {
        string id = p->identifier;
        string domain = CasProfileUtil::domain_table_fetch(p->identifier);
        fmt::print("{}\t{}\n", id, domain);
    }
}

void writeBinaryFile(const string &fileName, vector<uint64_t> arr) {
    ofstream f{fileName, ios::out};
    f.seekp(0);
    for (auto &k: arr)
        f << k << ' ';
    f.close();
}

vector<uint64_t> readBinaryFile(const string &fileName) {
    std::vector<uint64_t> ret;

    ifstream infile;
    infile.open(fileName, ios::in);
    std::string line;

    while (std::getline(infile, line)) {
        char *ptr; // declare a ptr pointer
        ptr = strtok((char *) line.c_str(), " "); // use strtok() function to separate string using comma (,) delimiter.
        // use while loop to check ptr is not null
        while (ptr != NULL) {
            ret.push_back(stoll(ptr, nullptr, 10));
            ptr = strtok(NULL, " ");
        }
    }

    infile.close();

    return ret;
}

vector<CasProfile *> CasProfileUtil::deserialize_profiles(std::filesystem::path directory) {
    vector<CasProfile *> profiles;
    for (const auto &entry: std::filesystem::directory_iterator(directory)) {
        if (entry.path().string().find("_mask") != string::npos) continue;

        CasProfile *profile = new CasProfile;

        string file_path = entry.path().string();
        phmap::BinaryInputArchive archive(file_path.c_str());
//        profile->hash_table.load(archive);
        profile->identifier = entry.path().stem().string();

        profile->binary_kmers = readBinaryFile(file_path);
        profile->binary_masks = readBinaryFile(file_path + "_mask");

        if (!CasProfileUtil::domain_table_contains(profile->identifier)) {
            fmt::print("deserialization of {} skipped due to domain table containment failure\n", profile->identifier);
            continue;
        }

        profiles.push_back(profile);
    }

    fmt::print("loaded {} profiles\n", profiles.size());
    return profiles;
}

CasProfile *profile_factory(const string &id, vector<string> sequences, ull k) {
    auto *profile = new CasProfile;
    profile->identifier = id;

    for (string &seq: sequences) {
        // clean seq
        seq.erase(remove(seq.begin(), seq.end(), '.'), seq.end());
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        for (ull i = 0; i < seq.size() - k + 1; i++) {
            string kmer = seq.substr(i, k);

            if (kmer.find('O') != string::npos || kmer.find('X') != string::npos || kmer.find('B') != string::npos ||
                kmer.find('Z') != string::npos)
                continue;

            try {
                auto b_kmer = Util::kmer_to_binary(kmer);

                profile->binary_kmers.push_back(b_kmer.k);
                profile->binary_masks.push_back(b_kmer.mask);

//                ui kmer_enc = Util::encode_amino_kmer(kmer);
//                profile->hash_table.insert(kmer_enc);
            }
            catch (const std::exception &e) {
                std::cerr << e.what() << '\n';
            }
        }
    }

    double sum = 0;
    for (ui i = 0; i < sequences.size(); i++) {
        sum += sequences[i].size();
        break; // Only consider first sequence
    }

    std::sort(sequences.begin(), sequences.end(),
              [](const std::string &first, const std::string &second) { return first.size() < second.size(); });
    profile->length_median = sequences[sequences.size() / 2].length();
    profile->length_mean = sum; // / sequences.size();
    profile->length_min = sequences[0].size();
    profile->length_max = sequences[sequences.size() - 1].size();

    fmt::print("{}\t{}\t{}\t{}\t{}\n", id, profile->length_median, profile->length_mean, profile->length_min,
               profile->length_max);

    return profile;
}

void lower_case_string(string &str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
}


// This needs to be versatile in that it can handle fasta seqs that are both single-line and multi-line
vector<CasProfile *> generate_from_fasta(std::filesystem::path fasta_dir) {
    vector<CasProfile *> profiles;
    vector<string> seqs;

    for (const auto &entry: std::filesystem::directory_iterator(fasta_dir)) {
        string identifier = entry.path().stem().string();

        if (!domain_map.contains(identifier)) {
            fmt::print("skipping generation of profile {} because the identifier is not in the domain map\n",
                       identifier);
            continue;
        }

        std::map<string, string> fasta_seqs = Util::parse_fasta(entry.path().string(), false);

        seqs.clear();
        for (auto const &seq: fasta_seqs) {
            string id = seq.first;
            boost::algorithm::to_lower(id);

            // We only want consensus sequence
            if (id.find("consensus") != string::npos)
                seqs.push_back(seq.second);
        }

        CasProfile *profile = profile_factory(identifier, seqs, CasProfileUtil::k);
        profiles.push_back(profile);
    }

    return profiles;
}

void CasProfileUtil::serialize(std::filesystem::path path_bin_pro) {
    auto serialize_profile = [&](CasProfile *profile) {
        std::filesystem::path file_name = path_bin_pro / profile->identifier;
//        phmap::BinaryOutputArchive archive(file_name.string().c_str());
//        profile->hash_table.dump(archive);

        writeBinaryFile(file_name, profile->binary_kmers);
        writeBinaryFile(file_name += "_mask", profile->binary_masks);
    };

    auto profiles = generate_from_fasta(Config::path_cas_profiles);
    std::for_each(profiles.begin(), profiles.end(), serialize_profile);
}
