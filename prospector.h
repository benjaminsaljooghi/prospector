typedef unsigned long long ull;
typedef unsigned int ui;
typedef signed long int li;
typedef unsigned char uc;

namespace Prospector
{
    const ui k_start = 20;
    const ui k_end = 60;
    const ui k_count = (Prospector::k_end - Prospector::k_start);
    const ui spacer_min =  21;
    const ui spacer_max = 72;
    const ui spacer_skip = (Prospector::spacer_min - 1);
    const ui repeats_min = 3;

    const ui size = 16;
    const ui map_size = 64;
    const ui bits = 2;
    const ui mutant_tolerance_ratio = 8;

    struct Encoding {
        ui* encoding;
        ui* d_encoding;
        ui size;
    };

    void device_init();
    Prospector::Encoding get_genome_encoding(const char* genome, ui genome_size);

    uc* get_qmap_small(ui* genome_encoding, ui genome_encoding_size);
    uc* get_qmap_big(const ui* genome_encoding, const ui genome_encoding_size, const ui* queries, const ui queries_size, ui map_size);
}

