// parser_cpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

string parse_fasta(string file_path)
{
    cout << "reading" << ++file_path << endl;

    string final_out = "";

    ifstream input(file_path);
    if (!input.good()) {
        cerr << "Error opening '" << file_path << "'. Bailing out." << endl;
        return "";
    }

    string line, name, content;
    while (getline(input, line)) {

        if (line.empty() || line[0] == '>') { // Identifier marker
            if (!name.empty()) { // Print out what we read from the last entry
                final_out += name + " : " + content + "\n";
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
        final_out += name + " : " + content + "\n";
    }

    return final_out;
}

int main(int argc, char** argv)
{


    if (argc <= 1) {
        cerr << "Usage: " << argv[0] << " [infile]" << endl;
        return -1;
    }

    string out = parse_fasta(argv[1]);
        
    cout << out;

    return 0;
}
