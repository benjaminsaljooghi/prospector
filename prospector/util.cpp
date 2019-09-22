#include "stdafx.h"
#include "util.h"

double Util::duration(clock_t begin)
{
	return (clock() - begin) / (double)CLOCKS_PER_SEC;
}

map<string, string> Util::parse_fasta(string file_path)
{
	cout << "reading: " << file_path << endl;
	ifstream input(file_path);
	if (!input.good())
	{
		throw runtime_error(strerror(errno));
	}

	map<string, string> seqs;
	string line, name, content;
	while (getline(input, line))
	{
		if (line.empty() || line[0] == '>') // Identifier marker
		{
			if (!name.empty())
			{
				// Get what we read from the last entry
				seqs[name] = content;
				name.clear();
			}
			if (!line.empty())
			{
				name = line.substr(1);
			}
			content.clear();
		}
		else if (!name.empty())
		{
			if (line.find(' ') != string::npos) // Invalid sequence--no spaces allowed
			{
				name.clear();
				content.clear();
			}
			if (line.find("\r") != string::npos)
			{
				throw runtime_error("File contains carriage returns (CRLF). Please reformat to LF.");
			}
			else
			{
				content += line;
			}
		}
	}
	if (!name.empty())
	{
		// Get what we read from the last 
		seqs[name] = content;
	}

	return seqs;
}

vector<int> Util::flatten(vector<vector<int>> vecs)
{
	vector<int> flattened;
	for (auto v : vecs)
		flattened.insert(flattened.end(), v.begin(), v.end());
	return flattened;
}

bool Util::subset(vector<int> a, vector<int> b)
{
	for (int e : a)
	{
		// Attempt to find e in in b. If the pointer arrives at b.end() then e is not inside b. Therefore. a is not contained within b
		if (find(b.begin(), b.end(), e) == b.end())
			return false;
	}
	// In all cases, elements within a were found in b. Therefore, a is a subset of b
	return true;
}


std::vector<std::string> Util::repeats(std::string genome, Util::Locus locus)
{
    std::vector<std::string> result;
    for (int index : locus.genome_indices)
    {
        result.push_back(genome.substr(index, locus.k).c_str());
    }
    return result;
}

std::vector<std::string> Util::spacers(std::string genome, Util::Locus locus)
{
    std::vector<std::string> result;
    for (unsigned int i = 0; i < locus.genome_indices.size()-1; i++)
    {
        int current_repeat_end = locus.genome_indices[i] + locus.k;
        int next_repeat_begin = locus.genome_indices[i+1];
        int spacer_size = next_repeat_begin - current_repeat_end;
        result.push_back(genome.substr(current_repeat_end, spacer_size));
    }   
    return result; 
}

std::string Util::seqs_to_fasta(vector<string> seqs)
{
    std::ostringstream string_stream;
    for (unsigned int i = 0; i < seqs.size(); i++)
    {
        string_stream << ">" << i << endl;
        string_stream << seqs[i] << endl;  
    }
    return string_stream.str();
}

