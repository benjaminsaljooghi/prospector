
#define DYAD_MIN 5

#define REPEAT_MIN 20
#define REPEAT_MAX 60

#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP 10

#define REPEATS_MIN 3
#define SCAN_DOMAIN 1000

#define ALLOW_DISCREPANT_LENGTHS false

#include <map>

const std::map<char, char> complements =
{
    { 'A', 'T' },
    { 'T', 'A' },
    { 'C', 'G' },
    { 'G', 'C' },
    { 'N', 'N' },
    { 'n', 'n' },
};

