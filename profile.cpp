
#include "stdafx.h"


class Profile
{
    public:
        string name;
        string seq;
        vector<string> kmers;
		
		Profile(string _name, string _path, unsigned int _k)
		{
			name = _name;
			seq = parse_fasta(_path).begin()->second;
			kmers = get_kmers(seq, _k);
		}
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
		
		
		ProfileExecution(Profile* _profile, Crispr* _crispr)
		{
			// how do we demarcate Cas genes?

			profile = _profile;
			crispr = _crispr;

			string a = "a";
			string b = "b";
			int a_b = a.compare(b);

			for (int query = 0; query < profile->kmers.size(); query++)
			{
				string query_kmer = profile->kmers[query];
				for (int target = 0; target < crispr->target_kmers.size(); target++)
				{
					// string target_kmer = crispr->target_kmers[target];
					// int comparison = query_kmer.compare(target_kmer);
					if (1 == 0)
					{
						locations_present[query_kmer].push_back(target);
						ordered_positions.push_back(target);
					}
				}
			}

			sort(ordered_positions.begin(), ordered_positions.end());

			hits = ordered_positions.size();
			hits_possible = (profile->kmers).size();
		}


		void print()
		{
			printf("profile %s; CRISPR %d %d", profile->name.c_str(), crispr->start, crispr->k);
			printf(": %zd/%zd\n", hits, hits_possible);
		}
};





