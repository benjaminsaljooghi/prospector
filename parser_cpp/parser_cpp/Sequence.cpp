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

CUDA_CALLABLE_MEMBER int Sequence::length()
{
    return seq.length();
}

CUDA_CALLABLE_MEMBER int Sequence::start()
{
    return start_pos;
}

CUDA_CALLABLE_MEMBER int Sequence::end()
{
    return start_pos + length() - 1;
}

CUDA_CALLABLE_MEMBER Sequence Sequence::subseq(int start, int length)
{
    return Sequence(seq.substr(start, length), start_pos + start);
}

CUDA_CALLABLE_MEMBER string Sequence::sequence()
{
    return seq;
}

CUDA_CALLABLE_MEMBER char Sequence::operator[](int i)
{
    return seq[i];
}

CUDA_CALLABLE_MEMBER bool Sequence::is_dyad()
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

CUDA_CALLABLE_MEMBER vector<Sequence> Sequence::dyads(int k)
{
    cout << "generating dyads for k: " << k << "... ";
    vector<Sequence> seqs;
    for (size_t i = 0; i < seq.length() - k + 1; i++)
    {
        Sequence kmer = subseq(i, k);
        if (kmer.is_dyad())
        {
            seqs.push_back(kmer);
        }
    }
    cout << "complete." << endl;
    return seqs;
}

CUDA_CALLABLE_MEMBER vector<Sequence> Sequence::dyads(int k_start, int k_end)
{
    vector<Sequence> seqs;
    for (int k = k_start; k < k_end; k++)
    {
        vector<Sequence> dyad_seqs = dyads(k);
        seqs.insert(seqs.end(), dyad_seqs.begin(), dyad_seqs.end());
    }
    return seqs;
}

CUDA_CALLABLE_MEMBER bool Sequence::operator<(const Sequence rhs) const
{
    bool comes_before = start_pos < rhs.start_pos;
    if (comes_before)
    {
        return true;
    }

    // sequences start at the same index, so compare based on size
    //bool smaller = seq.length() < rhs.seq.length();
    //if (smaller)
    //{
    //    return true;
    //}

    // sequences start at the same index and have the same size, so they must be equal
    return false;
}