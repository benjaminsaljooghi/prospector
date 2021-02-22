#pragma once

#include "stdafx.h"

#include "cas.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "array_discovery.h"

namespace Debug
{
	void visualize_map(string genome_path);

	void visualize_proximals(map<ull, vector<ull>> proximal_targets, string genome);
	
	vector<CasProfile*> cas_filter(vector<CasProfile*> profiles, string identifier);
	
	vector<Crispr*> crispr_filter(vector<Crispr*> crisprs, ull start, ull end);
	
	string translation_test(const string& genome, ull genome_start, ull genome_final, bool pos, ull debug_aminos);
	
	void translation_print(const string& genome, ull genome_start, ull genome_final, bool pos, ull debug_aminos);
	
	void triframe_print(const string& genome, ull genome_start, ull genome_final, bool pos);
	
	void cas_detect(const string& genome, ull genome_start, ull genome_final, bool pos, CasProfile* profile);
	
	void crispr_print(vector<Crispr*> crisprs, const string& genome, ull start, ull end);
	
	void genome_substr(const string& genome_path, ull genome_start, ull genome_final);

	void cartograph_interpreter(std::filesystem::path path, std::filesystem::path genome_dir);
}

