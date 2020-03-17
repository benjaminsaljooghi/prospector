#include "stdafx.h"


map <string, string> codon_table = {
    {"TTT", "F"},
    {"TTC", "F"},
    {"TTA", "L"},
    {"TTG", "L"},
    {"CTT", "L"},
    {"CTC", "L"},
    {"CTA", "L"},
    {"CTG", "L"},
    {"ATT", "I"},
    {"ATC", "I"},
    {"ATA", "I"},
    {"ATG", "M"},
    {"GTT", "V"},
    {"GTC", "V"},
    {"GTA", "V"},
    {"GTG", "V"},
    {"TCT", "S"},
    {"TCC", "S"},
    {"TCA", "S"},
    {"TCG", "S"},
    {"CCT", "P"},
    {"CCC", "P"},
    {"CCA", "P"},
    {"CCG", "P"},
    {"ACT", "T"},
    {"ACC", "T"},
    {"ACA", "T"},
    {"ACG", "T"},
    {"GCT", "A"},
    {"GCC", "A"},
    {"GCA", "A"},
    {"GCG", "A"},
    {"TAT", "Y"},
    {"TAC", "Y"},
    {"TAA", "STOP"},
    {"TAG", "STOP"},
    {"CAT", "H"},
    {"CAC", "H"},
    {"CAA", "Q"},
    {"CAG", "Q"},
    {"AAT", "N"},
    {"AAC", "N"},
    {"AAA", "K"},
    {"AAG", "K"},
    {"GAT", "D"},
    {"GAC", "D"},
    {"GAA", "E"},
    {"GAG", "E"},
    {"TGT", "C"},
    {"TGC", "C"},
    {"TGA", "STOP"},
    {"TGG", "W"},
    {"CGT", "R"},
    {"CGC", "R"},
    {"CGA", "R"},
    {"CGG", "R"},
    {"AGT", "S"},
    {"AGC", "S"},
    {"AGA", "R"},
    {"AGG", "R"},
    {"GGT", "G"},
    {"GGC", "G"},
    {"GGA", "G"},
    {"GGG", "G"}
};

map <char, char> complement_table = {
	{'A', 'T'},
	{'T', 'A'},
	{'C', 'G'},
	{'G', 'C'},
	{'N', 'N'},
	{'n', 'n'},
};



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
string reverse_complement (string seq)
{
    size_t len = seq.length();
    string complement = seq;

    for (size_t i = 0; i < len; i++)
    {
        complement[i] = complement_table[seq[i]];
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

		_count += upstream == complement_table[downstream] ? 0 : 1;
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
vector <string> get_kmers (string seq, unsigned int k)
{
	vector<string> kmers;
	for (size_t i = 0; i < seq.length() - k + 1; i++)
		kmers.push_back(seq.substr(i, k));
	return kmers;
}
map <string, string*> sixwaytranslation (string seq)
{

    size_t k = 3;

    auto amino_acid_seq = [&](string& dna, size_t start)
    {
        string __amino_acid_seq = "";
        for (size_t i = start; i + k < dna.size(); i += k)
        {
            string codon = dna.substr(i, k);
            string amino_acid = codon_table[codon];
            if (amino_acid != "STOP") __amino_acid_seq += amino_acid;
        }
        return __amino_acid_seq;
    };

    auto frame_shift = [&](string& domain)
    {
        vector<string> amino_acid_seqs = vector<string>(3);
        for (size_t frame = 0; frame < 3; frame++)
            amino_acid_seqs[frame] = amino_acid_seq(domain, frame);
        return &amino_acid_seqs[0];
    };

    string rc = reverse_complement(seq);
    string* pos_seqs = frame_shift(seq);
    string* neg_seqs = frame_shift(rc);

    map<string, string*> result = {
            {"pos", pos_seqs},
            {"neg", neg_seqs},
    };

    return result;
}
#undef LZZ_INLINE
