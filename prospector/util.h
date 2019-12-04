#include "stdafx.h"

namespace Util
{
	struct Locus
	{
		int k;
		std::vector<int> genome_indices;

		int start(int repeat)
		{
			return genome_indices[repeat];
		}
		int end(int repeat)
		{
			return genome_indices[repeat] + k;
		}

		bool repeat_substring(int _start, int _end, int repeat)
		{
			return _start >= start(repeat) && _end <= end(repeat);
		}

		bool repeat_substring(int _start, int _end)
		{
			for (int i = 0; i < genome_indices.size(); i++)
			{
				if (repeat_substring(_start, _end, i))
				{
					return true;
				}				
			}
			return false;
		}

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
