#pragma once
#include "stdafx.h"


namespace Util
{
	struct Locus
	{
		int k;
		std::vector<int> genome_indices;
	};

	struct Prospection
	{
		std::string genome;
		std::vector<Util::Locus> crisprs;
	};


	double duration(clock_t begin);
	map<string, string> parse_fasta(string file_path);
	vector<int> flatten(vector<vector<int>> vecs);
	bool subset(vector<int> a, vector<int> b);	
	vector<std::string> repeats(std::string, Util::Locus);
	vector<std::string> spacers(std::string, Util::Locus);
	std::string seqs_to_fasta(std::vector<std::string>);
	bool repeat_subset(Locus, Locus);
}
