
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
    ui genome_start;
    ui genome_final;
    bool pos;
    string raw;
    string pure;
    vector<string> pure_kmerized;
    vector<ui> pure_kmerized_encoded;
    vector<ll> pure_mapping;
};

struct Fragment
{
    string* reference_genome;

    Crispr* reference_crispr;
    Translation* reference_translation;
    CasProfile* reference_profile;

    vector<vector<ll>> clusters;
    ll clust_begin;
    ll clust_final;

    ll genome_begin;
    ll genome_final;

    ll expanded_genome_begin;
    ll expanded_genome_final;

    string to_string_debug()
    {
        string amino_domain = Util::translate_genome(*reference_genome, genome_begin, genome_final, reference_translation->pos);
        string amino_gene = Util::translate_genome(*reference_genome, expanded_genome_begin, expanded_genome_final, reference_translation->pos);

        string dna_domain = reference_genome->substr(genome_begin, genome_final - genome_begin);
        string dna_gene = reference_genome->substr(expanded_genome_begin, expanded_genome_final - expanded_genome_begin);

        string amino_buffer;

        ui begin_discrepant = (genome_begin - expanded_genome_begin);
        ui final_discpreant = (expanded_genome_final - genome_final);

        dna_gene.insert(begin_discrepant, "-");
        amino_gene.insert(begin_discrepant / 3, "-");

        dna_gene.insert(dna_gene.size() - final_discpreant, "-");
        amino_gene.insert(amino_gene.size() - (final_discpreant / 3), "-");

        std::ostringstream out;
        out << fmt::format("{}\t{}\n", reference_translation->pos ? "+" : "-", reference_profile->identifier);
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

    ui get_start()
    {
        return fragments[0]->expanded_genome_begin;
    }


    string to_string_debug()
    {
        return fmt::format("debugging info not implemented for MultiFragment struct");
    }

    string to_string_summary()
    {
        std::string result;
        
        ui min = fragments[0]->expanded_genome_begin;
        ui max = fragments[0]->expanded_genome_final;

        for (ui i = 1; i < fragments.size(); i++)
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
            result += fmt::format("{}-{} ", domain, id);
        }

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
    static const ui upstream_size = 20000;
    static const ui cluster_metric_min = 15;
    static const ui max_inter_cluster_dist = 2;

    vector<Translation*> get_triframe(const string& genome, ll genome_start, ll genome_final, bool pos);
    vector<Translation*> get_sixframe(const string& genome, ll genome_start, ll genome_final);
    vector<Translation*> crispr_proximal_translations(const string& genome, vector<Crispr*>&);
    vector<Fragment*> cas(vector<CasProfile*>& cas_profiles, vector<Translation*>&, string&);
}

