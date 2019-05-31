#include "Sequence.h"
#include "Crispr.h"

#include <string>

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
    std::sort(repeats.begin(), repeats.end(), Sequence::compare_start);
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