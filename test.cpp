
#include "stdafx.h"

set<string> get_bigrams(const string& str)
{
    set<string> bigrams;

    for (ull i = 0; i < str.size() - 1; i++)
    {
        bigrams.insert(str.substr(i, 2));
    }
    return bigrams;
}

double similarity(const string& a, const string& b)
{
    set<string> bigrams_a = get_bigrams(a);
    set<string> bigrams_b = get_bigrams(b);

    set<string> intersection;
    set_intersection(
        bigrams_a.begin(),
        bigrams_a.end(),
        bigrams_b.begin(),
        bigrams_b.end(),
        inserter(intersection, intersection.begin())
    );

    return (2.0 * intersection.size()) / (bigrams_a.size() + bigrams_b.size());
}


double spacer_similarity(const vector<string>& spacers)
{
    double total = 0;
    ull comps = 0;
    for (ull i = 0; i < spacers.size(); i++)
    {
        for (ull j = 0; j < i; j++)
        {
            total += similarity(spacers[i], spacers[j]);
            comps++;
        }
    }
    return total / (double) comps;
}

int main()
{

}

// # def get_bigrams(string):
// #     '''
// #     Takes a string and returns a list of bigrams
// #     '''
// #     s = string.lower()
// #     return {s[i:i+2] for i in range(len(s) - 1)}

// # def string_similarity(str1, str2):
// #     '''
// #     Perform bigram comparison between two strings
// #     and return a percentage match in decimal form
// #     '''
// #     pairs1 = get_bigrams(str1)
// #     pairs2 = get_bigrams(str2)
// #     print(f"{pairs1=} {pairs2=}")
// #     print(f"{pairs1 & pairs2}")
// #     return (2.0 * len(pairs1 & pairs2)) / (len(pairs1) + len(pairs2))

// # want_lower = [
// #     "CATCATAGGCGGAACTGGTAGGATGTACACAGCAAGTGGGTTTCC",
// #     "TTAATTGAATTATGAAATACAAGGTAATAACATATTTTGAGTTTCC",
// #     "TTTCTAGGAATGGGTAATTATAGCGAGCTAGAAAGCGTTTCC"
// # ]

// # want_higher = [
// #     "CATCATAGGCGGAACTGGTAGGATGTACACAGCAAGTGGGTTTCCG",
// #     "TTAATTGAATTATGAAATACAAGGTAATAACATATTTTGAGTTTCCG",
// #     "TTTCTAGGAATGGGTAATTATAGCGAGCTAGAAAGCGTTTCCG"
// # ]

// # def total_sim(strs):
// #     total = 0
// #     comps = 0
// #     for i in range(len(strs)):
// #         for j in range(i):
// #             print(f"{i=} {j=}")
// #             total += string_similarity(strs[i], strs[j])
// #             comps += 1
// #     return (float(total)) / (float(comps))


// # def total_sim(strs):
// #     total = 0
// #     comps = 0
// #     for i in range(len(strs)):
// #         for j in range(i):
// #             print(f"{i=} {j=}")
// #             total += string_similarity(strs[i], strs[j])
// #             comps += 1
// #     return (float(total)) / (float(comps))
        

// # out = total_sim(want_higher)

// # print(out)