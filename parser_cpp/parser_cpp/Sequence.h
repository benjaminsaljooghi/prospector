#pragma once

#include <vector>
#include <string>

class Sequence
{
    std::string seq;
    int start_pos;

    public:

    Sequence(std::string seq, int start);

    int length();

    int start();

    int end();

    Sequence subseq(int start, int length);

    std::string sequence();

    char operator[](int i);

    bool is_dyad();

    std::vector<Sequence> dyads(int k);

    std::vector<Sequence> dyads(int k_start, int k_end);

    bool operator<(const Sequence rhs) const;
};

