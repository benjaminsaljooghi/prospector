#include "Sequence.h"
#include "Crispr.h"

#include <string>

using namespace std;

Crispr::Crispr(int k)
{
    this->k = k;
}

void Crispr::add_repeat(int repeat)
{
    repeats.push_back(repeat);
}

int Crispr::last()
{
    return repeats.back();
}

void Crispr::sort_repeats()
{
    std::sort(repeats.begin(), repeats.end());
}

string Crispr::stringification(Sequence genome)
{
    sort_repeats();
    string str = genome.subseq(repeats[0], k).seq + ": ";
    for (int i = 0; i < repeats.size(); i++)
    {
        str += repeats[i] + " ";
    }
    return str;
}

bool Crispr::operator<(const Crispr& rhs) const
{
    int this_start = repeats.front();
    int rhs_start = rhs.repeats.front();
    bool index_before = this_start < rhs_start;
    if (index_before)
    {
        return true;
    }
    
    bool count_lower = repeats.size() < rhs.repeats.size();
    if (count_lower)
    {
        return true;
    }

    // crisprs start at the same position, and have the same repeat count. Therefore, let's compare based on sequences
    for (int i = 0; i < repeats.size(); i++)
    {
        int this_repeat = repeats[i];
        int other_repeat = rhs.repeats[i];
        if (this_repeat < other_repeat)
        {
            return true;
        }
    }

    if (k < rhs.k)
    {
        return true;
    }

    // this crispr is not < rhs
    return false;
}