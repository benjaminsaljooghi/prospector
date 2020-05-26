#pragma once


#include "stdafx.h"

#include "cas.h"
#include "crispr.h"
#include "util.h"

namespace Debug
{

	vector<CasProfile> cas_filter(vector<CasProfile> profiles, string gn);

	vector<Crispr> crispr_filter(vector<Crispr> crisprs, ui start, ui end);

	string translation_test(const string& genome, ui genome_start, ui genome_final, bool pos, ui debug_aminos);

	void crispr_print(vector<Crispr> crisprs, const string& genome, ui start, ui end);

}

