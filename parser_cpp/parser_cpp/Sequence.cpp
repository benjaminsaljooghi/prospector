#include "Sequence.h"
#include "consts.h"

#include <string>
#include <vector>
#include <map>

using namespace std;

Sequence::Sequence(string seq, int start)
{
    this->seq = seq;
    this->start_pos = start;
}

int Sequence::length()
{
    return seq.length();
}

int Sequence::start()
{
    return start_pos;
}

int Sequence::end()
{
    return start_pos + length() - 1;
}

Sequence Sequence::subseq(int start, int length)
{
    return Sequence(seq.substr(start, length), start_pos + start);
}

string Sequence::sequence()
{
    return seq;
}

char Sequence::operator[](int i)
{
    return seq[i];
}

bool Sequence::is_dyad()
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

vector<Sequence> Sequence::dyads(int k)
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

bool Sequence::compare_start(Sequence a, Sequence b)
{
    return a.start() < b.start();
}

