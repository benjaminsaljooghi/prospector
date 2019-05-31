#pragma once

#include "Sequence.h"

#include <vector>
#include <algorithm>

class Crispr
{
    public:

    std::vector<Sequence> repeats;

    void add_repeat(Sequence repeat);

    Sequence last();

    void sort_repeats();

    std::string to_string();
};