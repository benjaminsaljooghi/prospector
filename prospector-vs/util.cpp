#include "util.h"

string Util::translate_domain(const string& domain)
{
	ull codon_size = 3;
	string raw = "";
	for (ull i = 0; i + codon_size < domain.size(); i += codon_size)
	{
		string codon = domain.substr(i, codon_size);
		string amino_acid = codon_table.at(codon);
		raw += amino_acid;
	}
	return raw;
}

string Util::translate_genome(const string& genome, ui genome_start, ui genome_final, bool pos)
{
	string domain = genome.substr(genome_start, genome_final - genome_start);
	domain = pos ? domain : Util::reverse_complement(domain);
	return Util::translate_domain(domain);
}

bool Util::any_overlap(ui a_start, ui a_final, ui b_start, ui b_final)
{
	bool a_before_b = a_start <= b_start;
	bool b_before_a = b_start <= a_start;

	bool a_bleeds_into_b = a_before_b && a_final >= b_start;
	bool b_bleeds_into_a = b_before_a && b_final >= a_start;

	return a_bleeds_into_b || b_bleeds_into_a;
}


map<string, string> Util::parse_fasta (string file_path)
{
	// printf("reading %s... ", file_path.c_str());
	// auto start = time();
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

	// done(start);
	return seqs;
}


string Util::parse_fasta_single (string file_path)
{
	return Util::parse_fasta(file_path).begin()->second;
}


string Util::reverse_complement (string seq)
{
    ull len = seq.length();
    string complement = seq;

    for (ull i = 0; i < len; i++)
    {
        complement[i] = Util::complement_table.at(seq[i]);
    }

	reverse(complement.begin(), complement.end());
	return complement;
}
int Util::mismatch_count (string repeat)
{
	int _count = 0;

	ull k = repeat.size();
    unsigned int start_index = 0;
	unsigned int end_index = start_index + repeat.size() - 1;

	for (ull i = 0; i < k/2; i++)
	{
		char upstream = repeat[start_index + i];
		char downstream = repeat[end_index - i];

		_count += upstream == Util::complement_table.at(downstream) ? 0 : 1;
	}
	return _count;
}
string Util::seqs_to_fasta (vector <string> seqs)
{
    ostringstream string_stream;
    for (ull i = 0; i < seqs.size(); i++)
    {
        string_stream << ">" << i << endl;
        string_stream << seqs[i] << endl;  
    }
    return string_stream.str();
}

vector <string> Util::kmerize (string seq, ui k)
{
	vector<string> kmers;
	for (ui i = 0; i < seq.length() - k + 1; i++)
		kmers.push_back(seq.substr(i, k));
	return kmers;
}

ui Util::encode_amino_kmer(string kmer, ui k)
{
	ui encoded = 0;
	for (ui i = 0; i < k; i++)
		encoded += Util::amino_encoding.at(kmer[i]) << k * i;
	return encoded;
}

vector<ui> Util::encode_amino_kmers(vector<string> kmers, ui k)
{
	assert(kmers[0].size() == k);
	vector<ui> encoded(kmers.size());
	memset(&encoded[0], 0, sizeof(ui) * kmers.size());
	for (ui j = 0; j < kmers.size(); j++)
	{
		encoded[j] = encode_amino_kmer(kmers[j], k);

	}
	return encoded;
}



map<string, string> Util::load_genomes(string dir)
{
    map<string, string> genomes;
    for (const auto& entry : filesystem::directory_iterator(dir))
        genomes[entry.path().stem().string()] = (parse_fasta_single(entry.path().string()));
    return genomes;
}

