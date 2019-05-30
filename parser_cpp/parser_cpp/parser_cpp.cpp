using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <optional>
#include <functional>

#pragma region consts

#define DYAD_MIN 5

#define REPEAT_MIN 20
#define REPEAT_MAX 60

#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10

#define REPEATS_MIN 3
#define SCAN_DOMAIN 1000

#define ALLOW_DISCREPANT_LENGTHS false

const map<char, char> complements =
{
    { 'A', 'T' },
    { 'T', 'A' },
    { 'C', 'G' },
    { 'G', 'C' },
    { 'N', 'N' },
    { 'n', 'n' },
};

#pragma endregion


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

class Sequence
{
    string seq;
    int start_pos;

    public:

    Sequence(string seq, int start)
    {
        this->seq = seq;
        this->start_pos = start;
    }

    int length()
    {
        return seq.length();
    }

    int start()
    {
        return start_pos;
    }

    int end()
    {
        return start_pos + length() - 1;
    }

    Sequence subseq(int start, int length)
    {
        return Sequence(this->seq.substr(start, length), this->start + start);
    }

    char operator[](int i)
    {
        return seq[i];
    }

    bool is_dyad()
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

    vector<Sequence> dyads(int k)
    {
        vector<Sequence> seqs;
        for (size_t i = 0; i < seq.length() - k + 1; i++)
        {
            Sequence seq = seq.subseq(i, k);
            if (seq.is_dyad())
            {
                seqs.push_back(seq);
            }
        }
        return seqs;
    }


};

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

class Crispr
{
    public:

    vector<Sequence> repeats;

    void add_repeat(Sequence repeat)
    {
        repeats.push_back(repeat);
    }

};

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
        //try
        //{
            Sequence kmer = genome.subseq(index++, k);
            if (mutant(consensus, kmer))
            {
                crispr.add_repeat(kmer);
                index = kmer.start() + k + SPACER_SKIP;
                countdown = reset;
            }
        //}
        //catch (ArgumentOutOfRangeException)
        //{
            //cout << "index was out of bounds, continuing..." << endl;
            //break;
        //}

    }

    // Downstream scan
    index = consensus.start() - k - SPACER_SKIP;
    countdown = reset;
    while (countdown-- > 0)
    {
        //try
        //{
        Sequence kmer = genome.subseq(index--, k);
        if (mutant(consensus, kmer))
        {
            crispr.add_repeat(kmer);
            index = kmer.start() - k - SPACER_SKIP;
            countdown = reset;
        }
        //}
        //catch (ArgumentOutOfRangeException)
        //{
            //Console.WriteLine("Index was out of bounds. Continuing...");
            //break;
        //}
    }

    if (crispr.repeats.size() >= REPEATS_MIN)
    {
        return optional<Crispr>{crispr};
    }
    else
    {
        nullopt;
    }
}

vector<Crispr> discover_crisprs(Sequence genome, int k)
{
    vector<Crispr> crisprs;
    vector<Sequence> dyads = genome.dyads(k);


}


int main()
{

    string test_path = R"(P:\CRISPR\test_data\test.fasta)";
    string aureus_path = R"(P:\CRISPR\bacteria\aureus.fasta)";

    string aureus = parse_single_seq(aureus_path);
    cout << aureus.length() << endl;

    return 0;
}
