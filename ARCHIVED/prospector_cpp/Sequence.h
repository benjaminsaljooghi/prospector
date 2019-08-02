#pragma once

#include <vector>
#include <string>

class Sequence
{
    public:

    std::string seq;
    int start_pos;

    Sequence(std::string seq, int start);
    int length();
    int start();
    int end();
    Sequence subseq(int start, int length);
    char operator[](int i);
    bool dyad(int i, int k);
    bool dyad();
    std::vector<Sequence> dyads(int k);
    bool operator<(const Sequence rhs) const;
};

