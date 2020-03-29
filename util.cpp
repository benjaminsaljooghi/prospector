#include "util.h"


double duration (double begin)
{
  return omp_get_wtime() - begin;
}

void done (double start, string out)
{
	printf("%s done in %.3f seconds\n", out.c_str(), duration(start));
}

void done (double start)
{
	done(start, "");
}

bool subset (vector <int> a, vector <int> b)
{
	for (int e : a)
	{
		if (find(b.begin(), b.end(), e) == b.end())
			return false;
	}
	return true;
}



map <string, string> parse_fasta (string file_path)
{
	printf("reading %s... ", file_path.c_str());
	double start = omp_get_wtime();
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

	done(start);
	return seqs;
}


string parse_fasta_single (string file_path)
{
	return parse_fasta(file_path).begin()->second;
}


string reverse_complement (string seq)
{
    size_t len = seq.length();
    string complement = seq;

    for (size_t i = 0; i < len; i++)
    {
        complement[i] = complement_table.at(seq[i]);
    }

	reverse(complement.begin(), complement.end());
	return complement;
}
int mismatch_count (string repeat)
{
	int _count = 0;

	size_t k = repeat.size();
    unsigned int start_index = 0;
	unsigned int end_index = start_index + repeat.size() - 1;

	for (size_t i = 0; i < k/2; i++)
	{
		char upstream = repeat[start_index + i];
		char downstream = repeat[end_index - i];

		_count += upstream == complement_table.at(downstream) ? 0 : 1;
	}
	return _count;
}
string seqs_to_fasta (vector <string> seqs)
{
    ostringstream string_stream;
    for (size_t i = 0; i < seqs.size(); i++)
    {
        string_stream << ">" << i << endl;
        string_stream << seqs[i] << endl;  
    }
    return string_stream.str();
}

vector <string> kmerize (string seq, unsigned int k)
{
	vector<string> kmers;
	for (size_t i = 0; i < seq.length() - k + 1; i++)
		kmers.push_back(seq.substr(i, k));
	return kmers;
}




Translation::Translation(string& _seq, unsigned int k)
:
	nucleotide_sequence(_seq)
{
    size_t codon_size = 3;

    auto frame_shift = [&](string& dna)
    {
        vector<string> amino_acid_seqs{"", "", ""};
        for (size_t frame = 0; frame < 3; frame++)
		{
			for (size_t i = frame; i + codon_size < dna.size(); i += codon_size)
			{
				string codon = dna.substr(i, codon_size);
				string amino_acid = codon_table.at(codon);
				amino_acid_seqs[frame] += amino_acid;
			}
		}
		return amino_acid_seqs;
    };

    string rc = reverse_complement(nucleotide_sequence);
    vector<string> pos_seqs = frame_shift(nucleotide_sequence);
    vector<string> neg_seqs = frame_shift(rc);

    this->translations_raw["pos_0"] = pos_seqs[0];
	this->translations_raw["pos_1"] = pos_seqs[1];
	this->translations_raw["pos_2"] = pos_seqs[2];
	this->translations_raw["neg_0"] = neg_seqs[0];
	this->translations_raw["neg_1"] = neg_seqs[1];
	this->translations_raw["neg_2"] = neg_seqs[2];

	// the genome positions are cached and computed here, they are not computed on the fly
	// they are cached in the object via a kind of "map"
	// that is, if I ask for the "index" of a translation,
	// it gets first mapped to the index in teh translation_raw
	// and then that index is multiplied by 3 to get me 
	// the genome index

	for (auto const& [key, val] : this->translations_raw)
	{
		this->translations_pure[key] = "";
		
		size_t stop_count = 0;
		size_t index = 0;
		for (char elem : val)
		{
			if (elem == STOP_C)
			{
				stop_count++;
				continue;
			}
			this->translations_pure[key] += elem;
			pure_mapping[key].push_back(index + stop_count);
			index++;
		}


		this->translations_pure_kmerized[key] = kmerize(this->translations_pure[key], k);
	}


}


size_t Translation::frame_offset(string label)
{
	map<string, size_t> offset;

	offset["pos_0"] = 0;
	offset["pos_1"] = 1;
	offset["pos_2"] = 2;
	offset["neg_0"] = 0;
	offset["neg_1"] = 1;
	offset["neg_2"] = 2;

	return offset[label];

}