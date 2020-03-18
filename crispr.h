#pragma once

#include "stdafx.h"
#include "util.h"
#include "blast.h"


class Crispr
{
    public:

        // computed by constructor
        unsigned int* genome_indices;
        size_t size;
        unsigned int k;

        // computed by update
        unsigned int start;
        unsigned int end;
        vector<string> repeats;
        vector<string> spacers;
        double conservation_repeats;
        double conservation_spacers;
        double overall_heuristic;

        // computed by cache_upstream_kmers
        vector<string> target_kmers;

        Crispr(unsigned int, unsigned int*, unsigned int*);
        void update(string&);
        void print(string&, map<string, int>);
        bool operator>(const Crispr& obj);
        bool operator<(const Crispr&);
        void cache_upstream_kmers(string, size_t, unsigned int);
};

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

namespace CrisprUtil
{
    bool any_overlap(Crispr, Crispr);
    void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores);
    bool repeat_substring(Crispr b, unsigned int start, unsigned int end);
    bool repeat_subset(Crispr a, Crispr b);
    void cas(string genome, vector<Crispr> crisprs, const unsigned int k, const size_t upstream_size);
    vector<Crispr> get_domain_best(vector<Crispr> crisprs);
    vector<Crispr> score_filtered(vector<Crispr> crisprs, map<string, int> spacer_scores);
    void cache_crispr_information(vector<Crispr>& crisprs, string genome);
    void debug(string genome, vector<Crispr> crisprs);
    map<string, int> get_spacer_scores(vector<Crispr>& crisprs);
}
