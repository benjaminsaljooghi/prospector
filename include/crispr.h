#pragma once

#include "stdafx.h"
#include "util.h"
#include "locus.h"

class Crispr : public Locus
{
    public:

        // computed by constructor
        vector<ull> genome_indices;
        ull size;
        ull k;

        // computed by update
        ull genome_start;
        ull genome_final;
        vector<string> repeats;
        vector<string> spacers;
        double conservation_repeats;
        // double conservation_spacers;
        double conservation_spacers2;
        double spacer_variance;
        double overall_heuristic;

        // computed by cache_upstream_aminos
        // vector<string> target_kmers;

        Crispr(ull k, vector<ull> genome_indices, ull size);
        ~Crispr();
        void update(const string& genome);
        //void print_generic(const string& genome, function<void(string)>& print_spacer) const;
        //void print(const string& genome, map<string, int>) const;
        //void print(const string& genome) const;
        // void cache_upstream_aminos(string, ull, ull);
        string identifier_string() const;

        ull get_start();
        ull get_final();
        string to_string_debug();
        string to_string_summary();
        bool is_crispr() { return true; }
        bool is_gene() { return false; }
        bool is_domain() { return false; }
};


namespace CrisprUtil
{
    bool any_overlap(const Crispr* a, const Crispr* b);
    //void print(string genome, vector<Crispr*> crisprs, map<string, int> spacer_scores);
    //void print(string genome, vector<Crispr*> crisprs);
    bool repeat_substring(Crispr* b, ull start, ull end);
    bool repeat_subset(Crispr* a, Crispr* b);
    // void predict_cas(string genome, vector<Crispr> crisprs, const ull k, const ull upstream_size);
    
    vector<Crispr*> get_domain_best(vector<Crispr*> crisprs);
    vector<Crispr*> spacer_score_filtered(vector<Crispr*> crisprs, map<string, int> spacer_scores);
    
    void cache_crispr_information(const string& genome, vector<Crispr*>& crisprs);
    //map<string, int> get_spacer_scores(vector<Crispr*>& crisprs, string target_db_path);

    bool heuristic_less(const Crispr* a, const Crispr* b);
    bool heuristic_greater(const Crispr* a, const Crispr* b);

}

