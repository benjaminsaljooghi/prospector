#include "Sequence.h"
#include "consts.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

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

char Sequence::operator[](int i)
{
    return seq[i];
}

bool Sequence::dyad(int index, int k)
{
    for (int i = 0; i < DYAD_MIN; i++)
    {
        char beginning_upstream = seq[index];
        char end_downstream = seq[index + k - i - 1];
        char end_downstream_comp = complements.at(end_downstream);
        if (beginning_upstream != end_downstream_comp)
        {
            return false;
        }
    }
    return true;
}

bool Sequence::dyad()
{
    return dyad(0, length());
}

vector<Sequence> Sequence::dyads(int k)
{
    cout << "generating dyads for k: " << k << "... ";
    vector<Sequence> seqs;
    for (size_t i = 0; i < seq.length() - k + 1; i++)
    {
        Sequence kmer = subseq(i, k);
        if (kmer.dyad())
        {
            seqs.push_back(kmer);
        }
    }
    cout << "complete." << endl;
    return seqs;
}

bool Sequence::operator<(const Sequence rhs) const
{
    bool comes_before = start_pos < rhs.start_pos;
    if (comes_before)
    {
        return true;
    }

    // sequences start at the same index, so compare based on size
    bool smaller = seq.length() < rhs.seq.length();
    if (smaller)
    {
        return true;
    }

    // sequences start at the same index and have the same size, so they must be equal
    return false;
}