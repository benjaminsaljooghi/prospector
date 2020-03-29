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
        double spacer_variance;
        double overall_heuristic;

        // computed by cache_upstream_aminos
        // vector<string> target_kmers;

        Crispr(unsigned int, unsigned int*, unsigned int*);
        void update(string&);
        void print_generic(string& genome, function<void(string)>& print_spacer);
        void print(string&, map<string, int>);
        void print(string&);
        // void cache_upstream_aminos(string, size_t, unsigned int);
};

class CasProfile
{
    public:
        string name;
        vector<string> kmers;
		
		CasProfile(string _path, unsigned int _k);
};

class CrisprProfile
{
    public:

        const Crispr& crispr;
        const Translation& translation;
        CrisprProfile(Crispr& _crispr, Translation& _translation);
};



namespace CrisprUtil
{
    bool any_overlap(Crispr, Crispr);
    void print(string genome, vector<Crispr> crisprs, map<string, int> spacer_scores);
    void print(string genome, vector<Crispr> crisprs);
    bool repeat_substring(Crispr b, unsigned int start, unsigned int end);
    bool repeat_subset(Crispr a, Crispr b);
    // void cas(string genome, vector<Crispr> crisprs, const unsigned int k, const size_t upstream_size);
    vector<Crispr> get_domain_best(vector<Crispr> crisprs);
    vector<Crispr> spacer_score_filtered(vector<Crispr> crisprs, map<string, int> spacer_scores);
    void cache_crispr_information(vector<Crispr>& crisprs, string genome);
    void debug(string genome, vector<Crispr> crisprs);
    map<string, int> get_spacer_scores(vector<Crispr>& crisprs, string target_db_path);

    bool heuristic_less(const Crispr& a, const Crispr& b);
    bool heuristic_greater(const Crispr& a, const Crispr& b);
}

