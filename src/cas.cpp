#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include "cas_helpers.h"
#include "config.h"

vector<Fragment*> Cas::cas(vector<CasProfile*>& profiles, vector<Translation*>& translations, string& genome)
{
    vector<Fragment*> fragments;

    #pragma omp parallel for
    for (auto profile : profiles)
    {
        ui translation_index = -1;
        for (Translation* t : translations)
        {
            translation_index++;
            vector<ull> index;

            for (ull i = 0; i < t->pure_kmerized_encoded.size(); i++)
            {
                bool contains = profile->hash_table.contains(t->pure_kmerized_encoded[i]);
                if (contains)
                    index.push_back(i);
            }

            if (index.size() == 0)
                continue;

            // dump index to a file for analysis
            string file_name = Config::path_output / fmt::format("index_dump/{}_{}_{}_{}_{}_{}", CasProfileUtil::domain_table_fetch(profile->identifier), profile->identifier, t->genome_start, t->genome_final, translation_index, index.size());

            if (Config::dump_indices) {
                std::ofstream out(file_name);
                for (ull index_location : index)
                {
                    out << index_location << "\t" << t->genome_start + index_location << endl;
                }
                out.close();
            }

            vector<vector<ull>> clusters = cluster_index(index);
            if (!good_clusters(clusters))
            {
                continue;
            }
            else if (Config::dump_indices)
            {
                fmt::print("accepted {}\n", file_name);

                string file_name_clust = fmt::format("{}.clust", file_name);
                std::ofstream out(file_name_clust);
                for (vector<ull> clust : clusters)
                {
                    for (ull index_location : clust)
                    {
                        out << index_location << "\t" << t->genome_start + index_location << endl;
                    }
                    out << endl;
                }
                out.close();
            }


            Fragment* f = new Fragment;
            f->reference_genome = &genome;
            f->reference_crispr = t->reference_crispr;
            f->reference_translation = t;
            f->reference_profile = profile;
            f->clusters = clusters;

            compute_demarc(f);

            if (f->clust_final - f->clust_begin <= 15)
                continue;

            compute_details(f, genome);

            auto expansion = f->reference_translation->pos ? fragment_expansion_pos : fragment_expansion_neg;
            expansion(f, genome);

            if (f->genome_begin < 0)
            {
                std::exit(1);
            }

            // if (f->repetition_problem())
            // {
                // fmt::print("rejected on the basis of repetition problem\n");
                // continue;
            // }

            #pragma omp critical
            {
                fragments.push_back(f);
            }
        }
    }


    vector<Fragment*> filtered_fragments;

    for (Fragment* f : fragments)
    {
        if (CasProfileUtil::domain_table_contains(f->reference_profile->identifier))
            filtered_fragments.push_back(f);
        else
            delete f;
    }

    std::sort(filtered_fragments.begin(), filtered_fragments.end(), [](Fragment* a, Fragment* b) {return a->expanded_genome_begin < b->expanded_genome_begin; });

    return filtered_fragments;
}

vector<Translation*> Cas::get_triframe(const string& genome, ull genome_start, ull genome_final, bool pos)
{
    string domain = genome.substr(genome_start, genome_final - genome_start);
    //domain = pos ? domain : Util::reverse_complement(domain);

    if (!pos)
    {
        Util::reverse_complement(domain);
    }

    vector<Translation*> translations;
    for (ull frame = 0; frame < 3; frame++)
	{
        Translation* translation = new Translation;
        translation->pos = pos;
        translation->genome_start = pos ? genome_start + frame : genome_start;
        translation->genome_final = pos ? genome_final : genome_final - frame;

        translation->raw = Util::translate_domain(domain.substr(frame));

		translation->pure = "";

		ull stop_count = 0;
		ull index = 0;
		for (char elem : translation->raw )
		{
			if (elem == Util::stop_c)
			{
				stop_count++;
				continue;
			}
			translation->pure += elem;
			translation->pure_mapping.push_back(index + stop_count);
			index++;
		}
        translation->pure_kmerized = Util::kmerize(translation->pure, CasProfileUtil::k);
        translation->pure_kmerized_encoded = Util::encode_amino_kmers(translation->pure_kmerized, CasProfileUtil::k);
        translations.push_back(translation);
    }
    return translations;
}

vector <Translation*> Cas::get_sixframe(const string& genome, ull genome_start, ull genome_final)
{
    vector<Translation*> sixframe;

    for (Translation* translation : Cas::get_triframe(genome, genome_start, genome_final, true))
        sixframe.push_back(translation);

    for (Translation* translation : Cas::get_triframe(genome, genome_start, genome_final, false))
        sixframe.push_back(translation);

    return sixframe;
}

vector<Translation*> Cas::crispr_proximal_translations(const string& genome, vector<Crispr*>& crisprs)
{
    vector<Translation*> translations;


    std::sort(crisprs.begin(), crisprs.end(), [](Crispr* a, Crispr* b) { return a->genome_start < b->genome_start; });
    const ull min_translation_size = 20;
    const ull permitted_translation_overlap = 100;
    ull minimum_bound = 0;
    for (Crispr* c : crisprs)
    {
        ull translation_start;
        ull translation_final;

        // normal case
        translation_start = c->genome_start - Cas::upstream_size;
        translation_final = c->genome_final + Cas::upstream_size;

        // guard bounds
        translation_start = translation_start > c->genome_start ? 0 : translation_start; // underflow
        translation_final = std::min(translation_final, (ull)genome.size() - 1); // exceeds genome

        // bound the start
        translation_start = std::max(translation_start, minimum_bound);

        // require min size
        if (translation_final - translation_start < min_translation_size)
            continue;

        // new bound
        minimum_bound = translation_final - permitted_translation_overlap;
        assert(minimum_bound < translation_final); // guard underflow

        vector<Translation*> six_frame = Cas::get_sixframe(genome, translation_start, translation_final);
        for (Translation* t : six_frame) t->reference_crispr = c;
        translations.insert(translations.end(), six_frame.begin(), six_frame.end());
    }
    return translations;
}







#pragma clang diagnostic pop