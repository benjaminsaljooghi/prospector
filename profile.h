#include "stdafx.h"
#include "crispr.h"


class Profile
{
    public:
        string name;
        string seq;
        vector<string> kmers;
		
		Profile(string _name, string _path, unsigned int _k);
};

class ProfileExecution
{
    public:
        Profile* profile;
        Crispr* crispr;
        vector<int> ordered_positions;
        map<string, vector<int>> locations_present;
        size_t hits;
        size_t hits_possible;
		
		ProfileExecution(Profile* _profile, Crispr* _crispr);
		void print();
};





