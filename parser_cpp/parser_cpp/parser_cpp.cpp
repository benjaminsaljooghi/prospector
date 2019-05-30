#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

map<string,string> parse_fasta(string file_path)
{
    cout << "reading: " << file_path << endl;

    map<string, string> seqs;

    ifstream input(file_path);
    if (!input.good()) {
        throw "Error opening " + file_path;
    }

    string line, name, content;
    while (getline(input, line)) {

        if (line.empty() || line[0] == '>') { // Identifier marker
            if (!name.empty()) { // Print out what we read from the last entry
                //final_out += name + " : " + content + "\n";
                seqs[name] = content;
                name.clear();
            }
            if (!line.empty()) {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty()) {


            if (line.find(' ') != string::npos) { // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if (!name.empty()) { // Print out what we read from the last 
        seqs[name] = content;
    }

    return seqs;
}

string parse_single_seq(string file_path)
{
    map<string, string> seqs = parse_fasta(file_path);
    string seq = seqs.begin()->second;
    return seq;
}

int main()
{
    string test_path = R"(P:\CRISPR\test_data\test.fasta)";
    string aureus_path = R"(P:\CRISPR\bacteria\aureus.fasta)";
    string out = parse_single_seq(aureus_path);       
    cout << "the content is" << endl << out;
    return 0;
}
