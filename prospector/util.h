// C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// std
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

// My stuff
//#include "Sequence.h"

using namespace std;

// Timing

clock_t start;
#define BEGIN start = clock();
#define END printf("done in %.3f seconds\n", duration(start));

double duration(clock_t begin)
{
	return (clock() - begin) / (double)CLOCKS_PER_SEC;
}



map<string, string> parse_fasta(string file_path)
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


string parse_single_seq(string file_path)
{
	map<string, string> seqs = parse_fasta(file_path);
	string seq = seqs.begin()->second;
	//return Sequence(seq, 0);
	return seq;
}


string parse_genome(string file_path)
{
	printf("parse genome...\n");
	//string genome = parse_single_seq(file_path).seq;
	//return genome;
	return parse_single_seq(file_path);
}