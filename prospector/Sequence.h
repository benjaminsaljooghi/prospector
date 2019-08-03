#pragma once

#include <vector>
#include <string>

class Sequence
{
    public:

    std::string seq;
    int start_pos;

    Sequence();
    Sequence(std::string seq, int start);
    int length();
    int start();
    int end();
    Sequence subseq(int start, int length);
    char operator[](int i);
    bool is_dyad(int dyad_min);
    std::vector<Sequence> dyads(int dyad_min, int k);
    std::vector<Sequence> dyads(int dyad_min, int k_start, int k_end);
    bool operator<(const Sequence rhs) const;
};

