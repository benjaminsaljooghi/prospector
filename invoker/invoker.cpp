
// C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// C++
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include <ctime>
#include <string>

using namespace std;

vector<vector<int> > prospector_main(char** argv);

int main(int argc, char** argv)
{
    printf("running invoker\n");
    vector<vector<int> > crisprs = prospector_main(argv);

    for (auto vec : crisprs)
    {
        int k = vec[0];
        int genome_index = vec[1];
        printf("%d %d\n", genome_index, k);
    }
}