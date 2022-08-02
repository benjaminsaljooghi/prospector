//
// Created by zach on 27/07/22.
//

#include "config.h"
#include "util.h"
#include <getopt.h>
#include "fmt/format.h"

namespace Config
{
    const std::filesystem::path default_path_output = "./prospector_output";
    std::filesystem::path path_output;
    std::filesystem::path path_results;
    std::filesystem::path path_bin_pro;
    std::filesystem::path path_map_dom;
    std::filesystem::path path_genome;
    std::filesystem::path path_cas_profiles;

    int crispr_proximal_search = 0;
    int skip_serialisation = 0;
    int cas_only = 0;
    int dump_indices = 0;

    /**
     * Program Arguments
     * [REQUIRED] --prof=/path/to/cas_profiles/     : Path to directory containing Cas profiles for Cas detection
     * [REQUIRED] --domMap=/path/to/domain_map.tsv  : Path to domain map for Cas detection
     * [REQUIRED] --genome=/path/to/genome_folder/  : Path to folder containing genomes for analysis
     * [OPTIONAL] --out=/path/to/output_folder/     : Path to folder where all Prospector output will be stored, defaults to ./prospector_output/
     * [OPTIONAL] --skipSer                         : skip serialisation of profiles
     * [OPTIONAL] --proxSearch                      : Use the CRISPR Proximal Search method for finding Cas genes
     * [OPTIONAL] --casOnly                         : Only perform Cas detection, skipping CRISPR detection
     * [OPTIONAL] --dumpIndices                     : Dump each Cas detection index to file, in the output folder prospector_output/index_dump/
     */
    static struct option long_options[] = {
            {"prof", required_argument, nullptr,'p' },
            {"domMap", required_argument,nullptr,'d' },
            {"genome", required_argument,nullptr,'g' },
            {"out", required_argument,nullptr,'o' },
            {"skipSer", no_argument, &skip_serialisation,1},
            {"proxSearch", no_argument, &crispr_proximal_search,1 },
            {"casOnly", no_argument, &cas_only,1 },
            {"dumpIndices", no_argument, &dump_indices,1 },
            {0, 0, 0,0 },
    };
}


void print_help_message()
{
    fmt::print("-- prospector command line arguments --\n");
    fmt::print("[REQUIRED] --prof=/path/to/cas_profiles/     : Path to directory containing Cas profiles for Cas detection\n"
               "[REQUIRED] --domMap=/path/to/domain_map.tsv  : Path to domain map for Cas detection\n"
               "[REQUIRED] --genome=/path/to/genome_folder/  : Path to folder containing genomes for analysis\n"
               "[OPTIONAL] --out=/path/to/output_folder/     : Path to folder where all Prospector output will be stored, defaults to ./prospector_output/\n"
               "[OPTIONAL] --skipSer                         : skip serialisation of profiles\n"
               "[OPTIONAL] --proxSearch                      : Use the CRISPR Proximal Search method for finding Cas genes\n"
               "[OPTIONAL] --casOnly                         : Only perform Cas detection, skipping CRISPR detection\n"
               "[OPTIONAL] --dumpIndices                     : Dump each Cas detection index to file, in the output folder prospector_output/index_dump/\n");
}

void Config::parse_program_args(int argc, char *argv[])
{
    int opt;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv,"h", long_options, &long_index )) != -1) {
        switch (opt) {
            case 'h':
                print_help_message();
                exit(0);
            case 'p' :
                Config::path_cas_profiles = optarg;
                break;
            case 'd' :
                Config::path_map_dom = optarg;
                break;
            case 'g' :
                Config::path_genome = optarg;
                break;
            case 'o' :
                Config::path_output = optarg;
                break;
            default:
                break;
        }
    }

    Util::assert_file(Config::path_map_dom);
    Util::assert_file(Config::path_genome);
    Util::assert_file(Config::path_cas_profiles);

    // Output path is optional
    if (Config::path_output.empty()) Config::path_output = Config::default_path_output;
    Config::path_results = Config::path_output / "results";
    Config::path_bin_pro = Config::path_output / "serialised_profiles";

    // Create dirs in case they don't exist yet
    std::filesystem::create_directory(Config::path_output);
    std::filesystem::create_directory(Config::path_results);
    std::filesystem::create_directory(Config::path_bin_pro);
}
