using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <optional>
#include <functional>

#include "consts.h"
#include "Sequence.h"
#include "Crispr.h"


map<string, string> parse_fasta(string file_path)
{
    cout << "reading: " << file_path << endl;
    ifstream input(file_path);
    if (!input.good())
    {
        throw "Error opening " + file_path;
    }

    map<string, string> seqs;
    string line, name, content;
    while (getline(input, line))
    {
        if (line.empty() || line[0] == '>') // Identifier marker
        {
            if (!name.empty())
            {
                // Get what we read from the last entry
                seqs[name] = content;
                name.clear();
            }
            if (!line.empty())
            {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty())
        {
            if (line.find(' ') != string::npos) // Invalid sequence--no spaces allowed
            {
                name.clear();
                content.clear();
            }
            else
            {
                content += line;
            }
        }
    }
    if (!name.empty())
    {
        // Get what we read from the last 
        seqs[name] = content;
    }

    return seqs;
}

Sequence parse_single_seq(string file_path)
{
    map<string, string> seqs = parse_fasta(file_path);
    string seq = seqs.begin()->second;
    return Sequence(seq, 0);
}

vector<string> get_kmers(string sequence, int k)
{
    vector<string> seqs;
    for (size_t i = 0; i < sequence.length() - k + 1; i++)
    {
        seqs.push_back(sequence.substr(i, k));
    }
    return seqs;
}



bool mutant(Sequence a, Sequence b)
{
    if (!ALLOW_DISCREPANT_LENGTHS && a.length() != b.length())
    {
        throw exception();
    }

    int len = min(a.length(), b.length());

    int allowed_point_mutations = a.length() / 10;
    int point_mutations = 0;

    for (int i = 0; i < len; i++)
    {
        if (a[i] != b[i] && ++point_mutations > allowed_point_mutations)
        {
            return false;
        }
    }
    return true;
}

optional<Crispr> discover_crispr(Sequence genome, Sequence consensus)
{
    Crispr crispr;
    crispr.add_repeat(consensus);

    int k = consensus.length();

    // Upstream scan

    int index = consensus.start() + k + SPACER_SKIP;
    const int reset = SCAN_DOMAIN;
    int countdown = reset;
    while (countdown-- > 0)
    {
        Sequence kmer = genome.subseq(index++, k);
        if (mutant(consensus, kmer))
        {
            crispr.add_repeat(kmer);
            index = kmer.start() + k + SPACER_SKIP;
            countdown = reset;
        }
        if (index + k > genome.end())
        {
            break;
        }
    }

    // Downstream scan
    index = consensus.start() - k - SPACER_SKIP;
    countdown = reset;
    while (countdown-- > 0)
    {
        Sequence kmer = genome.subseq(index--, k);
        if (mutant(consensus, kmer))
        {
            crispr.add_repeat(kmer);
            index = kmer.start() - k - SPACER_SKIP;   
            countdown = reset;
        }
        if (index < genome.start())
        {
            break;
        }
    }

    if (crispr.repeats.size() >= 3)
    {
        cout << "CRISPR discovered" << endl;
        return optional<Crispr>{crispr};
    }
    else
    {
        cout << "nullopt discovered" << endl;
        return nullopt;
    }
}

vector<Crispr> discover_crisprs(Sequence genome, int k)
{
    cout << "discovering crisprs for k: " << k << endl;

    vector<Crispr> crisprs;
    vector<Sequence> dyads = genome.dyads(k);

    for (int i = 0; i < dyads.size(); i++)
    {
        Sequence dyad = dyads[i];
        optional<Crispr> crispr = discover_crispr(genome, dyad);
        if (crispr.has_value())
        {
            cout << "registering crispr..." << endl;
            crisprs.push_back(*crispr);
            i = (*crispr).last().end() + 1;
        }
    }
    return crisprs;
}

vector<Crispr> discover_crisprs(Sequence genome, int k_start, int k_end)
{

    vector<Crispr> all_crisprs;

    for (int k = k_start; k < k_end; k++)
    {
        vector<Crispr> crisprs = discover_crisprs(genome, k);
        all_crisprs.insert(all_crisprs.end(), crisprs.begin(), crisprs.end());
    }

    return all_crisprs;
}


int main()
{

    string test_path = R"(P:\CRISPR\test_data\test.fasta)";
    string aureus_path = R"(P:\CRISPR\bacteria\aureus.fasta)";
    string pyogenes_path = R"(P:\CRISPR\bacteria\pyogenes.fasta)";

    Sequence pyogenes = parse_single_seq(pyogenes_path);

    //vector<Crispr> crisprs = discover_crisprs(pyogenes, REPEAT_MIN, REPEAT_MAX);
    vector<Crispr> crisprs = discover_crisprs(pyogenes, 36);

    for (int i = 0; i < crisprs.size(); i++)
    {
        cout << crisprs[i].to_string() << endl;
    }

    return 0;
}
