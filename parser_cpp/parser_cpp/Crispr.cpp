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

string Crispr::to_string()
{
    sort_repeats();
    string str = "";
    for (int i = 0; i < repeats.size(); i++)
    {
        Sequence seq = repeats[i];
        str += seq.sequence() + " ";
    }
    return str;
}