#pragma once

#include "stdafx.h"

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

bool any_overlap(Crispr, Crispr);