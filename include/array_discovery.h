
#pragma once

#include "stdafx.h"
#include "prospector.h"
#include "util.h"
#include "crispr.h"
#include "debug.h"

namespace Array
{
	bool mutant(const char* genome, kmer* h, ui k, ull i, ull j, ui repeat_tolerance_ratio);
	vector<Crispr*> get_crisprs(string& genome);

}