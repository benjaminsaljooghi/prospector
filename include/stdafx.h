#include <cstdio>
#include <cstdlib>
#include <cstring>
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
#include <sstream>
#include <atomic>
#include <chrono>
#include <thread>
#include <cassert>
#include <omp.h>
//#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cmath>
#include <csignal>

#include <bitset>
#include <bit>
#include <cstddef>
#include <filesystem>

#include <chrono>

#define FMT_HEADER_ONLY 1

#include "fmt/core.h"
#include "fmt/format.h"
#include "time_util.h"


#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif

using namespace std;
namespace fs = std::filesystem;

typedef unsigned long long ull;
typedef long long ll;
typedef unsigned int ui;
typedef unsigned int kmer;