using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#define DYAD_MIN 5

#define REPEAT_MIN 20
#define REPEAT_MAX 60

#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10

#define REPEATS_MIN 3
#define SCAN_DOMAIN 1000



const map<char, char> complements =
{
    { 'A', 'T' },
    { 'T', 'A' },
    { 'C', 'G' },
    { 'G', 'C' },
    { 'N', 'N' },
    { 'n', 'n' },
};

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

string parse_single_seq(string file_path)
{
    map<string, string> seqs = parse_fasta(file_path);
    string seq = seqs.begin()->second;
    return seq;
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


bool is_dyad(string seq)
{
    int len = seq.length();
    for (int i = 0; i < DYAD_MIN; i++)
    {
        char beginning_upstream = seq[i];
        char end_downstream = seq[len - i - 1];
        char end_downstream_comp = complements.at(end_downstream);
        if (beginning_upstream != end_downstream_comp)
        {
            return false;
        }
    }
    return true;
}

vector<string> get_dyads(string sequence, int k)
{
    vector<string> seqs;
    for (size_t i = 0; i < sequence.length() - k + 1; i++)
    {
        string seq = sequence.substr(i, k);
        if (is_dyad(seq))
        {
            seqs.push_back(seq);
        }        
    }
    return seqs;
}

class Sequence
{
public:
    string seq;
    int start;

    Sequence(string seq, int start)
    {
        this->seq = seq;
        this->start = start;
    }

    int length()
    {
        return seq.length();
    }

    int start()
    {
        return start;
    }

    int end()
    {
        return start + length() - 1;
    }

    Sequence substring(int start, int length)
    {
        return Sequence(this->seq.substr(start, length), this->start + start);
    }
};

class Crispr
{
public:
    vector<Sequence> repeats;

};

Crispr discover_crispr(string genome, int consensus_index, int consensus_len)
{
    //Crispr crispr = new Crispr();
    //crispr.AddRepeat(consensus);
    string consensus = genome.substr(consensus_index, consensus_len))

    Cirpsr crispr;
    crispr.push_back(consensus)

    int k = consensus_size;
    int spacer_skip = 10;

    // Upstream scan
    int index = consensus.Start + k + spacer_skip;
    const int reset = SCAN_DOMAIN;
    int countdown = reset;
    while (countdown-- > 0)
    {
        try
        {
            Sequence kmer = genome.Substring(index++, k);
            if (Sequence.Mutant(consensus, kmer))
            {
                crispr.AddRepeat(kmer);
                index = kmer.Start + k + spacer_skip;
                countdown = reset;
            }
        }
        catch (ArgumentOutOfRangeException)
        {
            Console.WriteLine("Index was out of bounds. Continuing...");
            break;
        }

    }

    // Downstream scan
    index = consensus.Start - k - spacer_skip;
    countdown = reset;
    while (countdown-- > 0)
    {
        try
        {
            Sequence kmer = genome.Substring(index--, k);
            if (Sequence.Mutant(consensus, kmer))
            {
                crispr.AddRepeat(kmer);
                index = kmer.Start - k - spacer_skip;
                countdown = reset;
            }
        }
        catch (ArgumentOutOfRangeException)
        {
            Console.WriteLine("Index was out of bounds. Continuing...");
            break;
        }
    }

    return crispr.Repeats.Count >= REPEATS_MIN ? crispr : null;
}



}

int main()
{

    string test_path = R"(P:\CRISPR\test_data\test.fasta)";
    string aureus_path = R"(P:\CRISPR\bacteria\aureus.fasta)";

    string aureus = parse_single_seq(aureus_path);
    cout << aureus.length() << endl;

    vector<string> dyads = get_dyads(aureus, 6);
    cout << dyads.size() << endl;

    return 0;
}
