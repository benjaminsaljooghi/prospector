//
// Created by zach on 2/08/22.
//

#ifndef PROSPECTOR_SYSTEM_H
#define PROSPECTOR_SYSTEM_H

#include "locus.h"
#include "stdafx.h"

struct System
{
    vector<Locus*> loci;
    string type;

    ull crispr_count()
    {
        ull count = 0;
        for (Locus* l : loci)
        {
            if (l->is_crispr())
                count++;
        }
        return count;
    }

    ull cas_count()
    {
        ull count = 0;
        for (Locus* l : loci)
        {
            if (!l->is_crispr())
                count++;
        }
        return count;
    }

    void sort_loci()
    {
        sort(this->loci.begin(), this->loci.end(), [](Locus* a, Locus* b) { return a->get_start() - b->get_start(); });
    }

    ull get_start()
    {
        Locus* first = this->loci[0];
        return first->get_start();
    }

    ull get_final()
    {
        Locus* last = this->loci[this->loci.size()-1];
        return last->get_final();
    }

    string to_string_summary(string& genome_id)
    {
        std::ostringstream out;
        out << '>' << genome_id << endl;
        for (Locus* l : this->loci)
            out << l->to_string_summary() << endl;
        return out.str();
    }

    string to_string_debug(string& genome_id)
    {
        std::ostringstream out;
        out << '>' << genome_id << endl;
        for (Locus* l : this->loci)
            out << l->to_string_debug() << endl;
        return out.str();
    }

    bool legitimate_system()
    {
        return this->cas_count() >= 2;
    }

};

#endif //PROSPECTOR_SYSTEM_H
