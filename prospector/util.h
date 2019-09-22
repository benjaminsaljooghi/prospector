#pragma once
#include "stdafx.h"



#define PRINTF_BYTE_FORMAT_ALIGN 10

namespace Util
{
	double duration(clock_t begin);

	map<string, string> parse_fasta(string file_path);

	vector<int> flatten(vector<vector<int>> vecs);

	bool subset(vector<int> a, vector<int> b);

	struct Locus
	{
		int k;
		std::vector<int> genome_indices;
		// std::vector<int> spacer_scores;
	};

	struct Prospection
	{
		std::string genome;
		std::vector<Util::Locus> crisprs;
	};

	
	vector<std::string> repeats(std::string genome, Util::Locus locus);
	vector<std::string> spacers(std::string genome, Util::Locus locus);
	std::string seqs_to_fasta(std::vector<std::string> seqs);
}
