//
// Created by zach on 27/07/22.
//

#ifndef PROSPECTOR_CONFIG_H
#define PROSPECTOR_CONFIG_H

#include "path.h"

namespace Config
{
    extern const std::filesystem::path default_path_output;
    extern std::filesystem::path path_output;
    extern std::filesystem::path path_results;
    extern std::filesystem::path path_bin_pro;

    extern std::filesystem::path path_map_dom;
    extern std::filesystem::path path_genome;
    extern std::filesystem::path path_cas_profiles;

    extern int crispr_proximal_search;
    extern int skip_serialisation;
    extern int dump_indices;
    extern int cas_threshold;
    extern int cas_chunk_length;

    void parse_program_args(int argc, char *argv[]);
}


#endif //PROSPECTOR_CONFIG_H
