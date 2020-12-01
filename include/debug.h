#pragma once

#include "stdafx.h"

#include "cas.h"
#include "crispr.h"
#include "util.h"
#include "prospector.h"
#include "array_discovery.h"

namespace Debug
{
	void visualize_map(string& genome_path);

	void visualize_proximals(map<ui, vector<ui>> proximal_targets, string genome);
	
	vector<CasProfile*> cas_filter(vector<CasProfile*> profiles, string identifier);
	
	vector<Crispr*> crispr_filter(vector<Crispr*> crisprs, ui start, ui end);
	
	string translation_test(const string& genome, ui genome_start, ui genome_final, bool pos, ui debug_aminos);
	
	void translation_print(const string& genome, ui genome_start, ui genome_final, bool pos, ui debug_aminos);
	
	void triframe_print(const string& genome, ui genome_start, ui genome_final, bool pos);
	
	void cas_detect(const string& genome, ui genome_start, ui genome_final, bool pos, CasProfile* profile);
	
	void crispr_print(vector<Crispr*> crisprs, const string& genome, ui start, ui end);
	
	void genome_substr(const string& genome_path, ui genome_start, ui genome_final);

	void sage_interpreter(string path, string genome_dir);
}

