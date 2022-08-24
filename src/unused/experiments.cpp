//
// Created by zach on 2/08/22.
//

// // hardcoded against cas1 for now
// vector<ui> hammer(vector<Fragment*> fragments) {

//     // runtime:
//     // - write candidate sequences to file
//     // - system call of hmmscan
//     // - parse results
//     // - make decisions based on results

//     vector<string> seqs;
//     for (Fragment* fragment : fragments) {
//         seqs.push_back(fragment->protein_sequence);
//     }

//     string fasta = Util::seqs_to_fasta(seqs);

//     // write candidate sequence to file
//     std::ofstream out("/home/ben/crispr/hmm/candidates_cas1.txt");
//     out << fasta;
//     out.close();

//     // system call of hmmscan
//     string call = "hmmscan --tblout tbloutty.txt /home/ben/crispr/hmm/Cas_Cas1.hmm /home/ben/crispr/hmm/candidates_cas1.txt > /dev/null";
//     system(call.c_str());

//     // parse tblout
//     std::ifstream infile("tbloutty.txt");
//     string line;
//     vector<ui> acceptable_indices;

//     while (std::getline(infile, line))
//     {
//         // fmt::print("line: {}\n", line);
//         if (line.starts_with("#"))
//             continue;

//         vector<string> parse_result = Util::parse(line, " ");
//         string query_name = parse_result[2];
//         double score = std::stod(parse_result[5]);
//         fmt::print("{} {}\n", query_name, score);
//         acceptable_indices.push_back(std::stoi(query_name));
//     }

//     for (ui ting : acceptable_indices)
//     {
//         fmt::print("{}\n", ting);
//     }

//     return acceptable_indices;
// }

// vector<Fragment*> acceptable_indices_filter_methodology(vector<Fragment*> raw_fragments)
// {
//     // perform an analysis here of the fragments (either here or multifragments, not sure yet, but let's focus on proof of concept first)
//     // let's focus just on cas1 for now
//     vector<ui> acceptable_indices = hammer(raw_fragments);

//     fmt::print("fasta enumeration of raw fragments:\n");
//     for (ui i = 0; i < raw_fragments.size(); i++) {
//         fmt::print("{}\n", i);
//     }

//     fmt::print("we now have acceptable indices\n");
//     for (ui index : acceptable_indices) {
//         fmt::print("{}\n", index);
//     }

//     vector<Fragment*> fragments;
//     for (ui i = 0; i < raw_fragments.size(); i++) {
//         if (Util::contains(acceptable_indices, i))
//             fragments.push_back(raw_fragments[i]);
//     }

//     return fragments;
// }



//map<string, System*> system_map;
//for (Fragment* f : fragments)
//{
//    if (!f->is_gene())
//    {
//        continue;
//    }

//    string id = f->reference_crispr->identifier_string();
//    if (system_map.contains(id))
//    {
//        system_map[id]->fragments.push_back(f);
//        continue;
//    }

//    System* system = new System;
//    system->crispr = f->reference_crispr;
//    system->fragments.push_back(f);
//    system_map[id] = system;
//}

//vector<System*> systems;

//for (map<string, System*>::iterator i = system_map.begin(); i != system_map.end(); i++)
//{
//    systems.push_back(i->second);
//}

//for (System* s : systems)
//{
//    s->type = resolve_type(s->fragments);
//}

//std::sort(systems.begin(), systems.end(), [](System* a, System* b) {return a->get_start() - b->get_start(); });
//for (System* s : systems)
//out_gene << s->to_string();


//std::map<string, std::set<string>> load_type_table(std::filesystem::path path)
//{
//    std::map<string, std::set<string>> type_table;
//
//    std::ifstream table(path);
//
//    if (!table)
//    {
//        throw runtime_error(path.string());
//    }
//
//    string line;
//    while (getline(table, line))
//    {
//        auto split = Util::parse(line, "\t");
//        auto type = split[0];
//        auto required = split[1];
//        auto additional_any = split[2];
//
//        std::set<string> genes;
//
//        for (string gene : Util::parse(, ","))
//            genes.insert(gene);
//
//        type_table[links] = genes;
//    }
//}

//string resolve_type(vector<Fragment*> fragments)
//{
//    return "N";
//    for (Fragment* f : fragments)
//    {
//        string r = CasProfileUtil::domain_table_fetch(f->reference_profile->identifier);
//        for (const auto& [a, b] : type_table)
//            if (b.contains(r))
//                return a;
//    }
//    return "N";
//}

//vector<string> get_all_hmm_files()
//{
//    vector<string> paths;
//
//    // step 1: iterate over the dir:
//    std::filesystem::path directory = "/home/ben/crispr/prospector-data/bin_hmm_makarova";
//
//    // does the split[0] of the stem match the cas_type in a case insensitive way?
//    for (const auto& entry : std::filesystem::directory_iterator(directory))
//    {
//        string filename = entry.path().filename().string();
//        vector<string> extensions = Util::parse(filename, ".");
//        string last_extension = extensions[extensions.size()-1];
//
//        if (last_extension != "hmm")
//            continue;
//
//        string cas_type_of_file = Util::parse(filename, "_")[0];
//        paths.push_back(entry.path().string());
//    }
//    return paths;
//}
//
//vector<string> get_hmm_files(string cas_type)
//{
//    vector<string> paths;
//
//    // step 1: iterate over the dir:
//    std::filesystem::path directory = "/home/ben/crispr/prospector-data/bin_hmm_makarova";
//
//    // does the split[0] of the stem match the cas_type in a case insensitive way?
//    for (const auto& entry : std::filesystem::directory_iterator(directory))
//    {
//        string filename = entry.path().filename().string();
//        vector<string> extensions = Util::parse(filename, ".");
//        string last_extension = extensions[extensions.size()-1];
//
//        if (last_extension != "hmm")
//            continue;
//
//        string cas_type_of_file = Util::parse(filename, "_")[0];
//        if (boost::iequals(cas_type_of_file, cas_type))
//        {
//            paths.push_back(entry.path().string());
//        }
//    }
//    return paths;
//}
//
//string get_hmm_file(string hmm_identifier)
//{
//    fmt::print("seeking HMM file with identifier: {}\n", hmm_identifier);
//    return fmt::format("/home/ben/crispr/hmm/hmm_files/{}.HMM", hmm_identifier);
//}
//
//bool have_hmm_file(string hmm_file)
//{
//    bool exists = std::filesystem::exists(hmm_file);
//    if (!exists)
//    {
//        fmt::print("does NOT exist: {}\n", hmm_file);
//    }
//    return exists;
//}

//bool singleton_hammer(string& singleton_seq, string hmm_file)
//{
//    if (!have_hmm_file(hmm_file))
//    {
//        throw new exception();
//    }
//
//
//    vector<string> singleton_vec;
//    singleton_vec.push_back(singleton_seq);
//    string fasta = Util::seqs_to_fasta(singleton_vec);
//
//    std::filesystem::path temp_fasta_file = "/media/ben/temp/crispr_lfs/temp_fasta.txt";
//
//    // write candidate sequence to file
//    std::ofstream out(temp_fasta_file.string());
//    out << fasta;
//    out.close();
//
//    string call = fmt::format("hmmscan --cpu 10 --tblout tbloutty.txt {} {} > /dev/null", hmm_file, temp_fasta_file.string());
//    system(call.c_str());
//
//    // parse tblout
//    std::ifstream infile("tbloutty.txt");
//    string line;
//    vector<ui> acceptable_indices;
//
//    while (std::getline(infile, line))
//    {
//        // fmt::print("line: {}\n", line);
//        if (line.starts_with("#"))
//            continue;
//
//        vector<string> parse_result = Util::parse(line, " ");
//        string query_name = parse_result[2];
//        double score = std::stod(parse_result[5]);
//        fmt::print("{} {}\n", query_name, score);
//        acceptable_indices.push_back(std::stoi(query_name));
//    }
//
//    // return true;
//    assert(acceptable_indices.size() <= 1);
//    return acceptable_indices.size() > 0;
//}

// bool singleton_hammer(Fragment* fragment) {
//     string hmm_file;

//     // attempt an exact match
//     hmm_file = get_hmm_file(fragment->reference_profile->identifier);
//     if (!have_hmm_file(hmm_file)) {
//         // exact match could be found, so this is probably like a COG profile. In which case we get the profile's domain and then we find a suitable HMM for that domain.
//     }

//     return singleton_hammer(fragment->protein_sequence, hmm_file);
// }

//bool multimatch_hammer(Fragment* fragment) {
//    string cas_type_of_fragment = CasProfileUtil::domain_table_fetch(fragment->reference_profile->identifier);
//    vector<string> hmm_files = get_hmm_files(cas_type_of_fragment);
//
//    for (string hmm_file : hmm_files) {
//        if (singleton_hammer(fragment->protein_sequence, hmm_file)) {
//            return true;
//        }
//    }
//    return false;
//}

// vector<Fragment*> individuated_singleton_methodology(vector<Fragment*> raw_fragments) {
//     vector<Fragment*> fragments;
//     for (Fragment* fragment : raw_fragments) {
//         if (singleton_hammer(fragment)) {
//             fragments.push_back(fragment);
//         }
//         else {
//             fmt::print("rejecting {} {} {} {} on the basis of hammer\n", fragment->reference_profile->identifier, CasProfileUtil::domain_table_fetch(fragment->reference_profile->identifier), fragment->expanded_genome_begin, fragment->expanded_genome_final);
//         }
//     }
//     return fragments;
// }

//vector<Fragment*> multimatch_singleton_methodology(vector<Fragment*> raw_fragments) {
//    vector<Fragment*> fragments;
//    for (Fragment* fragment : raw_fragments) {
//        if (multimatch_hammer(fragment)) {
//            fragments.push_back(fragment);
//        }
//        else {
//            fmt::print("rejecting {} {} {} {} on the basis of hammer\n", fragment->reference_profile->identifier, CasProfileUtil::domain_table_fetch(fragment->reference_profile->identifier), fragment->expanded_genome_begin, fragment->expanded_genome_final);
//        }
//    }
//    return fragments;
//}

// void hmm_db_experiment()
// {
//     vector<string> hmm_files = get_all_hmm_files();

//     std::filesystem::path path_prodigal_proteins = "/home/ben/crispr/prospector-util/my.proteins.faa";

//     map<string, string> proteins = Util::parse_fasta(path_prodigal_proteins, false);

//     map<string, vector<Hit*>> protein_id_to_hits;

//     for (auto & [identifier, sequence] : proteins)
//     {
// 		sequence.erase(remove(sequence.begin(), sequence.end(), '*'), sequence.end());
//     }

//     for (auto & [identifier, sequence] : proteins)
//     {
//         fmt::print("{}\n", identifier);
//         for (string& hmm_file : hmm_files)F
//         {
//             bool positive = singleton_hammer(sequence, hmm_file);
//             if (positive) {
//                 fmt::print("got a positive on {}\n", sequence);
//             }
//         }
//     }

// }

//struct Hit
//{
//    CasProfile* profile;
//    ui hit_count;
//};

//bool claim_validation(string& protein_sequence, string cas_claim)
//{
//    vector<string> hmm_files = get_hmm_files(cas_claim);
//    fmt::print("executing {} hmm files\n", hmm_files.size());
//    for (string hmm_file : hmm_files)
//    {
//        if (singleton_hammer(protein_sequence, hmm_file))
//        {
//            return true;
//        }
//    }
//    return false;
//}
//
//
//string hit_assessment(Hit* hit, string& protein_sequence)
//{
//    string cas_claim = CasProfileUtil::domain_table_fetch(hit->profile->identifier);
//    bool validated = claim_validation(protein_sequence, cas_claim);
//    return fmt::format("{} {} {} {}\n", hit->profile->identifier, cas_claim, hit->hit_count, validated);
//
//}
//
//void analyze_prodigal_proteins(vector<CasProfile*>& profiles)
//{
//    auto start_prodigal = time();
//
//    std::filesystem::path path_prodigal_proteins = "/home/ben/crispr/prospector-util/my.proteins.faa";
//
//    map<string, string> proteins = Util::parse_fasta(path_prodigal_proteins, false);
//
//    map<string, vector<Hit*>> protein_id_to_hits;
//
//    for (auto & [identifier, sequence] : proteins)
//    {
//        sequence.erase(remove(sequence.begin(), sequence.end(), '*'), sequence.end());
//
//        auto kmers = Util::kmerize(sequence, 6);
//        auto kmers_encoded = Util::encode_amino_kmers(kmers, 6);
//
//        vector<Hit*> hits;
//
//        for (signed long p = 0; p < profiles.size(); p++)
//        {
//            CasProfile* profile = profiles[p];
//            // vector<ull> index;
//            ull count = 0;
//            for (ull i = 0; i < kmers_encoded.size(); i++)
//            {
//                bool contains = profile->hash_table.contains(kmers_encoded[i]);
//                if (contains)
//                    count++;
//                // index.push_back(i);
//            }
//
//            Hit* hit = new Hit;
//            hit->profile = profile;
//            // hit->hit_count = index.size();
//            hit->hit_count = count;
//
//            hits.push_back(hit);
//        }
//
//        // sort hits
//        Util::sort(hits, [](Hit* a, Hit* b) { return a->hit_count > b->hit_count; });
//
//        // store hits
//        protein_id_to_hits[identifier] = hits;
//    }
//
//    start_prodigal = time(start_prodigal, "hit counts");
//
//    std::ofstream outtie("outtie.txt");
//
//    for (auto const& [identifier, hits] : protein_id_to_hits)
//    {
//        string protein_sequence = proteins[identifier];
//        outtie << fmt::format("{}\n", identifier);
//        outtie << hit_assessment(hits[0], protein_sequence);
//        outtie << hit_assessment(hits[1], protein_sequence);
//        outtie << hit_assessment(hits[2], protein_sequence);
//        outtie << "\n\n";
//    }
//
//    start_prodigal = time(start_prodigal, "hmmer");
//
//    outtie.close();
//
//    start_prodigal = time(start_prodigal, "full prodigal analysis");
//}


//Gene* gene_from_fragments(vector<Fragment*>& fragments)
//{
//    sort(fragments.begin(), fragments.end(), [](const Fragment* a, const Fragment* b) { return a->genome_begin < b->genome_begin;  });
//
//    Gene* g = new Gene;
//    g->gn = fragments[0]->reference_profile->gn;
//    g->fragments = fragments;
//    return g;
//}

//vector<Gene> best_genes(vector<Gene> genes)
//{
//    map<string, Gene> types;
//    for (const Gene& gene : genes)
//    {
//        auto gn = gene.fragments[0].reference_profile->gn;
//        if (types.contains(gn) && gene.size() <= types.at(gn).size())
//            continue;
//        types[gn] = gene;
//    }
//    vector<Gene> __genes;
//    for (auto pairing : types)
//        __genes.push_back(pairing.second);
//
//    return __genes;
//}

//vector<Gene*> genes_from_fragments(vector<Fragment*>& fragments)
//{
//    map<string, vector<Fragment*>> gene_fragments;
//    for (Fragment* a : fragments)
//        gene_fragments[a->reference_profile->gn].push_back(a);
//
//    vector<Gene*> genes;
//    for (auto gene : gene_fragments)
//    {
//        genes.push_back(gene_from_fragments(gene.second));
//    }
//
//    //genes = best_genes(genes);
//    sort(genes.begin(), genes.end(), [](Gene* a, Gene* b) { return a->fragments[0]->genome_begin < b->fragments[0]->genome_begin; } );
//
//    return genes;
//}


//map<string, vector<Gene*>> Cas::assemble_genes(const vector<Crispr*>& crisprs, const vector<Fragment*>& fragments)
//{
//    //map<string, Crispr> crispr_map;
//
//    map<string, vector<Fragment*>> crispr_fragments;
//    map<string, vector<Gene*>> crispr_genes;
//
//
//    for (const Crispr* c : crisprs)
//    {
//        string c_string = c->identifier_string();
//        for (Fragment* f : fragments)
//        {
//            if (f->reference_crispr->start == c->start && f->reference_crispr->k == c->k)
//            {
//                crispr_fragments[c_string].push_back(f);
//            }
//        }
//    }
//
//    for (auto [c_string, fragments] : crispr_fragments)
//    {
//        crispr_genes[c_string] = genes_from_fragments(fragments);
//    }
//
//    return crispr_genes;
//}

