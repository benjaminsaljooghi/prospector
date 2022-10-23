//
// Created by zach on 27/07/22.
//

#include "config.h"
#include "util.h"
#include <getopt.h>
#include "fmt/format.h"

#define PROF_OPTARG 'p'
#define DOM_MAP_OPTARG 'd'
#define GENOME_OPTARG 'g'
#define OUT_OPTARG 'o'
#define CAS_THRESHOLD_OPTARG 't'
#define CAS_CHUNK_LENGTH_OPTARG 'l'

namespace Config {
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
    int cas_threshold = 50; // Default to 50
    int cas_chunk_length = 1000; // Default to 1000

    // Program Arguments
    static struct option long_options[] = {
            {"prof",           required_argument, nullptr, PROF_OPTARG},
            {"domMap",         required_argument, nullptr, DOM_MAP_OPTARG},
            {"genome",         required_argument, nullptr, GENOME_OPTARG},
            {"out",            required_argument, nullptr, OUT_OPTARG},
            {"casThreshold",   required_argument, nullptr, CAS_THRESHOLD_OPTARG},
            {"casChunkLength", required_argument, nullptr, CAS_CHUNK_LENGTH_OPTARG},
            {"skipSer",        no_argument,       &skip_serialisation,     1},
            {"proxSearch",     no_argument,       &crispr_proximal_search, 1},
            {"casOnly",        no_argument,       &cas_only,               1},
            {"dumpIndices",    no_argument,       &dump_indices,           1},
            {0, 0,                                0,                       0},
    };
}


void print_help_message() {
    fmt::print("-- prospector command line arguments --\n");
    fmt::print(
            "\t[REQUIRED] --prof=/path/to/cas_profiles/     : Path to directory containing Cas profiles for Cas detection\n"
            "\t[REQUIRED] --domMap=/path/to/domain_map.tsv  : Path to domain map for Cas detection\n"
            "\t[REQUIRED] --genome=/path/to/genome_folder/  : Path to folder containing genomes for analysis\n"
            "\t[OPTIONAL] --out=/path/to/output_folder/     : Path to folder where all Prospector output will be stored, defaults to ./prospector_output/\n"
            "\t[OPTIONAL] --skipSer                         : skip serialisation of profiles\n"
            "\t[OPTIONAL] --proxSearch                      : Use the CRISPR Proximal Search method for finding Cas genes\n"
            "\t[OPTIONAL] --casOnly                         : Only perform Cas detection, skipping CRISPR detection\n"
            "\t[OPTIONAL] --dumpIndices                     : Dump each Cas detection index to file, in the output folder prospector_output/index_dump/\n"
            "\t[OPTIONAL] --casThreshold=50                 : Threshold score for a given Cas gene prediction to be considered\n"
            "\t[OPTIONAL] --casChunkLength=1000             : Length (in bp) of genome chunk to be analysed for Cas gene prediction at a time\n");
}

void Config::parse_program_args(int argc, char *argv[]) {
    int opt;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'h':
                print_help_message();
                exit(0);
            case PROF_OPTARG :
                Config::path_cas_profiles = optarg;
                break;
            case DOM_MAP_OPTARG :
                Config::path_map_dom = optarg;
                break;
            case GENOME_OPTARG :
                Config::path_genome = optarg;
                break;
            case OUT_OPTARG :
                Config::path_output = optarg;
                break;
            case CAS_THRESHOLD_OPTARG :
                Config::cas_threshold = stoi(optarg);
                break;
            case CAS_CHUNK_LENGTH_OPTARG :
                Config::cas_chunk_length = stoi(optarg);
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
