
#include "util.h"

std::vector<std::string> Util::parse(std::string str, std::string delim)
{
	//std::string s = "scott>=tiger>=mushroom";
	//std::string delimiter = ">=";

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
	//std::cout << s << std::endl;
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

string Util::translate_genome(const string& genome, ui genome_begin, ui genome_final, bool pos)
{
	string domain = genome.substr(genome_begin, genome_final - genome_begin);
	if (!pos)
	{
		Util::reverse_complement(domain);
	}
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


string Util::load_genome(string path)
{

	fmt::print("reading {}\n", path);


	ifstream input(path);
	if (!input.good())
	{
		throw runtime_error(strerror(errno));
	}


	auto check_line = [](string& line) {
		assert(line.find(' ') == string::npos);
		assert(line.find("\r") == string::npos);
	};


	string line;
	string name;
	string content;

	getline(input, line);
	assert(line.starts_with('>'));
	name = line.substr(1);
	
	while (true)
	{
		getline(input, line);

		if (line.empty() || line[0] == '>')
		{

			content.erase(std::remove(content.begin(), content.end(), 'N'), content.end());
			return content;

		}

		check_line(line);
		content += line;
	}

	assert(false);
}

//map<string, string> Util::parse_fasta (string file_path)
//{
//	// printf("reading %s... ", file_path.c_str());
//	// auto start = time();
//	ifstream input(file_path);
//	if (!input.good())
//	{
//		throw runtime_error(strerror(errno));
//	}
//
//	map<string, string> seqs;
//	string line, name, content;
//	while (getline(input, line))
//	{
//		if (line.empty() || line[0] == '>') // Identifier marker
//		{
//			if (!name.empty())
//			{
//				// Get what we read from the last entry
//				seqs[name] = content;
//				name.clear();
//			}
//			if (!line.empty())
//			{
//				name = line.substr(1);
//			}
//			content.clear();
//		}
//		else if (!name.empty())
//		{
//			if (line.find(' ') != string::npos) // Invalid sequence--no spaces allowed
//			{
//				name.clear();
//				content.clear();
//			}
//			if (line.find("\r") != string::npos)
//			{
//				throw runtime_error("File contains carriage returns (CRLF). Please reformat to LF.");
//			}
//			else
//			{
//				content += line;
//			}
//		}
//	}
//	if (!name.empty())
//	{
//		// Get what we read from the last 
//		seqs[name] = content;
//	}
//
//	// done(start);
//	return seqs;
//}

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

		_count += upstream == Util::complement(downstream) ? 0 : 1;
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

ui Util::encode_amino_kmer(const string& kmer)
{
	ull k = kmer.size();
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
		encoded[j] = encode_amino_kmer(kmers[j]);
	return encoded;
}



//map<string, string> Util::load_genomes(string dir)
//{
//	auto start = time();
//    map<string, string> genomes;
//    for (const auto& entry : filesystem::directory_iterator(dir))
//        genomes[entry.path().stem().string()] = (load_genome(entry.path().string()));
//	start = time(start, "genome load");
//	return genomes;
//}
//
