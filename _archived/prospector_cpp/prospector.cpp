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
#include <set>
#include <ctime>


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



//bool mutant(Sequence a, Sequence b)
//{
//    if (!ALLOW_DISCREPANT_LENGTHS && a.length() != b.length())
//    {
//        throw exception();
//    }
//
//    int len = min(a.length(), b.length());
//
//    int allowed_point_mutations = a.length() / 10;
//    int point_mutations = 0;
//
//    for (int i = 0; i < len; i++)
//    {
//        if (a[i] != b[i] && ++point_mutations > allowed_point_mutations)
//        {
//            return false;
//        }
//    }
//    return true;
//}

bool mutant(int a, int b, int k, string* genome)
{
    int allowed_point_mutations = k / 10;
    int point_mutations = 0;
    for (int i = 0; i < k; i++)
    {
        if ((*genome)[a + i] != (*genome)[b + i] && ++point_mutations > allowed_point_mutations)
        {
            return false;
        }
    }
    return true;
}

optional<Crispr> discover_crispr(string genome, int i, int k, int start, int end)
{

    Crispr crispr = Crispr(k);
    crispr.add_repeat(i);

    // Upstream scan
    int index = i + k + SPACER_SKIP;
    const int reset = SCAN_DOMAIN;
    int countdown = reset;

   
    while (countdown-- > 0)
    {
        if (index + k > end)
        {
            break;
        }
        if (mutant(i, index, k, &genome))
        {
            crispr.add_repeat(index);
            index += k + SPACER_SKIP;
            countdown = reset;
        }
        else
        {
            index++;
        }
    }

    // Downstream scan
    index = i + k + SPACER_SKIP;
    countdown = SCAN_DOMAIN;
    while (countdown-- > 0)
    {
        if (index < start)
        {
            break;
        }
        if (mutant(i, index, k, &genome))
        {
            crispr.add_repeat(index);
            index -= k + SPACER_SKIP;
            countdown = SCAN_DOMAIN;
        }
        else
        {
            index--;
        }
    }

    if (crispr.repeats.size() >= REPEATS_MIN)
    {
        crispr.sort_repeats();
        return optional<Crispr>{crispr};
    }

    return nullopt;
}

set<Crispr> discover_crisprs(Sequence genome, int k)
{
    cout << "discovering crisprs for k: " << k << endl;
    set<Crispr> crisprs;
    //vector<Sequence> dyads = genome.dyads(k);
    int start = genome.start();
    int end = genome.end();
    string seq = genome.seq;
    for (int i = 0; i < genome.length() - k + 1; i++)
    {
        //cout << i;
        //Sequence dyad = dyads[i];
        //cout << "\rexamining dyad " << i << "/" << dyads.size()-1 << " with start " << dyad.start() << "/" << genome.length();
        
        if (!genome.dyad(i, k))
            continue;
        
        optional<Crispr> crispr = discover_crispr(seq, i, k, start, end);
        if (crispr.has_value())
        {
            //cout << " -> CRISPR discovered at consensus start " << dyad.start() << endl;
            crisprs.insert(*crispr);
        }
    }
    cout << endl;
    return crisprs;
}

set<Crispr> discover_crisprs(Sequence genome, int k_start, int k_end)
{
    set<Crispr> all_crisprs;
    for (int k = k_start; k < k_end; k++)
    {
        set<Crispr> crisprs = discover_crisprs(genome, k);
        all_crisprs.insert(crisprs.begin(), crisprs.end());
    }
    return all_crisprs;
}


int main()
{


    //Crispr a = Crispr(35);
    //Crispr b = Crispr(35);


    //a.add_repeat(1684293);
    //a.add_repeat(1684293);
    //a.add_repeat(1684419);

    //b.add_repeat(1684245);
    //b.add_repeat(1684245);
    //b.add_repeat(1684371);
    //b.add_repeat(1684497);


    //cout << (a < b) << endl;
    //cout << (b < a) << endl;


    //return 0;



    clock_t start;
    double duration;
    start = clock();

    string test_path = R"(P:\CRISPR\data\test.fasta)";
    string aureus_path = R"(P:\CRISPR\data\aureus.fasta)";
    string pyogenes_path = R"(P:\CRISPR\data\pyogenes.fasta)";
    Sequence pyogenes = parse_single_seq(pyogenes_path);
    set<Crispr> crisprs = discover_crisprs(pyogenes, 35, 36);

    cout << "discovered CRISPRs: " << endl;
    for (auto c : crisprs)
    {
        cout << c.stringification(pyogenes) << endl;
    }

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "printf: " << duration << '\n';

    return 0;
}
