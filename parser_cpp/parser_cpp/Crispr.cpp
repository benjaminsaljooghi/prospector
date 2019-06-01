#include "Sequence.h"
#include "Crispr.h"

#include <string>
#include <algorithm>

using namespace std;

void Crispr::add_repeat(Sequence repeat)
{
    repeats.push_back(repeat);
}

Sequence Crispr::last()
{
    return repeats.back();
}

void Crispr::sort_repeats()
{
    std::sort(repeats.begin(), repeats.end());
}

string Crispr::stringification()
{
    sort_repeats();
    string str = repeats[0].sequence() + ": ";
    for (int i = 0; i < repeats.size(); i++)
    {
        int repeat_start = repeats[i].start();
        string repeat_str = to_string(repeat_start);
        str += repeat_str + " ";
    }
    return str;
}

bool Crispr::operator<(const Crispr& rhs) const
{

    Sequence this_start = repeats.front();
    Sequence rhs_start = rhs.repeats.front();
    //bool index_before = this_start.start() < rhs_start.start();
    int this_start_pos = this_start.start();
    int rhs_start_pos = rhs_start.start();

    if (this_start_pos < rhs_start_pos)
    {
        return true;
    }
    
    bool count_lower = repeats.size() < rhs.repeats.size();
    if (count_lower)
    {
        return true;
    }

    // crisprs start at the same position, and have the same repeat count. Therefore, let's compare based on sequences
    for (int i = 0; i < rhs.repeats.size(); i++)
    {
        Sequence this_repeat = repeats[i];
        Sequence other_repeat = rhs.repeats[i];
        if (this_repeat < other_repeat)
        {
            return true;
        }
    }

    // this crispr is not < rhs
    return false;

}