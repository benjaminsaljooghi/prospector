
#pragma once

#include "stdafx.h"
#include "prospector.h"
#include "util.h"
#include "crispr.h"
#include "debug.h"

namespace Array
{
	bool mutant(const char* genome, ui* h, ui k, ui i, ui j, ui repeat_tolerance_ratio);
	vector<Crispr*> get_crisprs(string& genome);

}