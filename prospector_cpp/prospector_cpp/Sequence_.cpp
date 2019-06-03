#include "Sequence.h"


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

    string sequence()
    {
        return seq;
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