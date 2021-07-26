
#pragma once

#include "stdafx.h"
#include <filesystem>
#include "fmt/core.h"
#include "fmt/format.h"


#include "util.h"
#include "crispr.h"
#include "cas_profiles.h"
#include "debug.h"
#include "locus.h"


namespace fs = std::filesystem;

struct Translation
{
    Crispr* reference_crispr;
    ull genome_start;
    ull genome_final;
    bool pos;
    string raw;
    string pure;
    vector<string> pure_kmerized;
    vector<kmer> pure_kmerized_encoded;
    vector<ull> pure_mapping;
};

struct Fragment
{
    string* reference_genome;

    Crispr* reference_crispr;
    Translation* reference_translation;
    CasProfile* reference_profile;

    string protein_sequence;

    vector<vector<ull>> clusters;
    ull clust_begin;
    ull clust_final;

    ull genome_begin;
    ull genome_final;

    ull expanded_genome_begin;
    ull expanded_genome_final;

    

    bool repetition_problem()
    {
        auto kmers = Util::kmerize(this->protein_sequence, CasProfileUtil::k); // no need to do this sicnce it's already been kmerized etc but keep it for now
        const int domain = 15;
        const int contiguous_check_requirement = 2;
        const int last_checkable_query = kmers.size() - domain;
        for (int i = 0; i < last_checkable_query; i++)
        {
            int how_many_repetitions = 0;
            for (int j = 0; j < domain; j++)
            {
                if (kmers[i] == kmers[j])
                    how_many_repetitions++;
            }
            if (how_many_repetitions >= contiguous_check_requirement)
                return true;
        }
        return false;
    } 

    string to_string_debug()
    {
        string amino_domain = Util::translate_genome(*reference_genome, genome_begin, genome_final, reference_translation->pos);
        string amino_gene = Util::translate_genome(*reference_genome, expanded_genome_begin, expanded_genome_final, reference_translation->pos);

        string dna_domain = reference_genome->substr(genome_begin, genome_final - genome_begin);
        string dna_gene = reference_genome->substr(expanded_genome_begin, expanded_genome_final - expanded_genome_begin);

        string amino_buffer;

        ull begin_discrepant = (genome_begin - expanded_genome_begin);
        ull final_discpreant = (expanded_genome_final - genome_final);

        dna_gene.insert(begin_discrepant, "-");
        amino_gene.insert(begin_discrepant / 3, "-");

        dna_gene.insert(dna_gene.size() - final_discpreant, "-");
        amino_gene.insert(amino_gene.size() - (final_discpreant / 3), "-");

        std::ostringstream out;
        out << fmt::format("{}\t{}\t{}\n", reference_translation->pos ? "+" : "-", reference_profile->identifier, CasProfileUtil::domain_table_fetch(reference_profile->identifier));
        out << fmt::format("\t{}...{}\n", genome_begin, genome_final);
        out << fmt::format("\t{}...{}\n", expanded_genome_begin, expanded_genome_final);
        out << fmt::format("\t{}\n", amino_gene);
        out << fmt::format("\t{}\n", dna_gene);
        return out.str();
    }

};

struct MultiFragment : public Locus
{
    vector<Fragment*> fragments;

    virtual ~MultiFragment()
    {
        for (Fragment* f : this->fragments)
        {
            delete f;
        }
    }


    ull get_start()
    {
        return fragments[0]->expanded_genome_begin;
    }

    ull get_final()
    {
        return fragments[fragments.size()-1]->expanded_genome_final;
    }


    string to_string_debug()
    {
        std::ostringstream out;
        out << "-------------------------" << endl;
        for (Fragment* f : fragments)
        {
            out << f->to_string_debug();
            out << endl;
        }
        out << "-------------------------" << endl;
        return out.str();
    }

    string to_string_summary()
    {
        std::string result;
        
        ull min = fragments[0]->expanded_genome_begin;
        ull max = fragments[0]->expanded_genome_final;

        for (ull i = 1; i < fragments.size(); i++)
        {
            if (fragments[i]->expanded_genome_begin < min)
                min = fragments[i]->expanded_genome_begin;

            if (fragments[i]->expanded_genome_final > max)
                max = fragments[i]->expanded_genome_final;
        }

        result = fmt::format("{}\t{}\t{}\t", 
            //fragments[0]->expanded_genome_begin,
            //fragments[0]->expanded_genome_final,
            min,
            max,
            fragments[0]->reference_translation->pos ? "+" : "-"
        );

        for (Fragment* fragment : fragments)
        {
            string id = fragment->reference_profile->identifier;
            string domain = CasProfileUtil::domain_table_contains(id) ? CasProfileUtil::domain_table_fetch(id) : id;
            result += fmt::format("{},", domain);
        }

        result = result.substr(0, result.size() - 1);

        result += "\t";

        for (Fragment* fragment : fragments)
        {
            string id = fragment->reference_profile->identifier;
            result += fmt::format("{},", id);
        }

        result = result.substr(0, result.size() - 1);

        return result;
    }

    bool is_crispr()
    {
        return false;
    }

    bool is_domain()
    {
        return !CasProfileUtil::domain_table_contains(fragments[0]->reference_profile->identifier);
    }

    bool is_gene()
    {
        return !is_domain();
    }
};



namespace Cas
{
    static const ull upstream_size = 30000;
    static const ull cluster_metric_min = 5; // lower = more sensitive
    static const ull max_inter_cluster_dist = 15; // higher = more sensitive

    vector<Translation*> get_triframe(const string& genome, ull genome_start, ull genome_final, bool pos);
    vector<Translation*> get_sixframe(const string& genome, ull genome_start, ull genome_final);
    vector<Translation*> crispr_proximal_translations(const string& genome, vector<Crispr*>&);
    vector<Fragment*> cas(vector<CasProfile*>& cas_profiles, vector<Translation*>&, string&);
}

