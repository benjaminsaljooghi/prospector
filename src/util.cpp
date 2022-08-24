
#include "util.h"

std::vector<std::string> Util::parse(std::string str, std::string delim)
{
	std::vector < std::string> tokens;

	size_t pos = 0;
	std::string token;
	while ((pos = str.find(delim)) != std::string::npos) {

		if (pos == 0) {
			pos = str.find_first_not_of(delim);
			str.erase(0, pos);
			continue;
		}

		token = str.substr(0, pos);
		//std::cout << token << std::endl;
		str.erase(0, pos + delim.length());

		if (token != delim)
			tokens.push_back(token);
	}
	tokens.push_back(str);
	return tokens;
}

std::vector<std::string> Util::parse_exact(std::string str, std::string delim)
{
	std::vector < std::string> tokens;

	size_t pos = 0;
	std::string token;
	while ((pos = str.find(delim)) != std::string::npos) {

		token = str.substr(0, pos);
		//std::cout << token << std::endl;
		str.erase(0, pos + delim.length());

		tokens.push_back(token);
	}
	tokens.push_back(str);
	return tokens;
}


ui Util::difference_cpu(const ui& _a, const ui& _b)
{
	ui _xor = (_a ^ _b);
	ui evenBits = _xor & 0xAAAAAAAAAAAAAAAAull;
	ui oddBits = _xor & 0x5555555555555555ull;
	ui comp = (evenBits >> 1) | oddBits;
	return __builtin_popcount(comp);
}

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

string Util::translate_genome(const string& genome, ull genome_begin, ull genome_final, bool pos)
{
	string domain = genome.substr(genome_begin, genome_final - genome_begin);
	if (!pos)
	{
		Util::reverse_complement(domain);
	}
	return Util::translate_domain(domain);
}

bool Util::any_overlap(ull a_start, ull a_final, ull b_start, ull b_final)
{
	bool a_before_b = a_start <= b_start;
	bool b_before_a = b_start <= a_start;

	bool a_bleeds_into_b = a_before_b && a_final >= b_start;
	bool b_bleeds_into_a = b_before_a && b_final >= a_start;

	return a_bleeds_into_b || b_bleeds_into_a;
}


Util::GenomeIdSequenceMap Util::load_genome(const std::filesystem::path& path)
{
	fmt::print("Reading {}...\n", path.string());

	ifstream input(path);
	if (!input.good())
	{
		throw runtime_error("Input file is not good!");
	}

    return parse_fasta(path.string(), true);
}

void remove_rna(string& content)
{
    content.erase(std::remove_if(content.begin(), content.end(), [](char c) { return !(c == 'A' || c == 'C' || c == 'G' || c == 'T'); }), content.end());
}

// dna = true for dna
// dna = false for amino
Util::GenomeIdSequenceMap Util::parse_fasta(const string& file_path, bool dna)
{
	ifstream input(file_path);
	if (!input.good())
	{
		throw runtime_error("input not good!");
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
				if (dna) remove_rna(content);
                seqs[name] = content;

				name.clear();
			}

			if (line[0] == '>')
			{
                // Genome name formatted like '>GENOME_ID <genome_sequence>'.
                // Find first occurrence of ' ' and remove it by subtracting 1.
                // Start search at position 1 to remove '>' from the name.
				name = line.substr(1, line.find(' ') - 1);

                if (!seqs[name].empty()) fmt::print("Duplicate genome found in file: {}\n", name);
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
		if (dna) remove_rna(content);
		seqs[name] = content;
	}

	return seqs;
}

char Util::complement(char nuc)
{
	switch (nuc)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	default:
		return 'N';
	}
}

void Util::complement(string& seq)
{
	for (char& c : seq) c = Util::complement(c);
	
}

void Util::reverse_complement(string& seq)
{
	complement(seq);
	reverse(seq.begin(), seq.end());
}


ull Util::mismatch_count (string repeat)
{
	ull _count = 0;

	ull k = repeat.size();
    ull start_index = 0;
	ull end_index = start_index + repeat.size() - 1;

	for (ull i = 0; i < k/2; i++)
	{
		char upstream = repeat[start_index + i];
		char downstream = repeat[end_index - i];

		_count += upstream == Util::complement(downstream) ? 0 : 1;
	}
	return _count;
}
string Util::seqs_to_fasta(vector<string>& seqs)
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
	for (ull i = 0; i < seq.length() - k + 1; i++)
		kmers.push_back(seq.substr(i, k));
	return kmers;
}

kmer Util::encode_amino_kmer(const string& seq)
{
	ull k = seq.size();
	kmer encoded = 0;
	for (ui i = 0; i < k; i++)
		encoded += Util::amino_encoding.at(seq[i]) << k * i;
	return encoded;
}

vector<kmer> Util::encode_amino_kmers(vector<string> kmers, ull k)
{
	assert(kmers[0].size() == k);
	vector<kmer> encoded(kmers.size());
	memset(&encoded[0], 0, sizeof(kmer) * kmers.size());
	for (ull j = 0; j < kmers.size(); j++)
		encoded[j] = encode_amino_kmer(kmers[j]);
	return encoded;
}

void Util::assert_file(std::filesystem::path path)
{
    if (!std::filesystem::exists(path))
    {
        throw std::filesystem::filesystem_error("Could not locate filepath: ", path, error_code());
    }
}
