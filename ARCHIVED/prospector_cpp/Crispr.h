#pragma once

#include "Sequence.h"

#include <vector>
#include <algorithm>

class Crispr
{
    public:

    int k;
    std::vector<int> repeats;

    Crispr(int k);
    void add_repeat(int repeat);
    int last();
    void sort_repeats();
    std::string stringification(Sequence genome);
    bool operator<(const Crispr& rhs) const;
};