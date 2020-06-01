typedef unsigned long long ull;
typedef unsigned int ui;
typedef signed long int li;
typedef unsigned char uc;

namespace Prospector
{
    static const ui k_start = 20;
    static const ui k_end = 60;
    static const ui k_count = (Prospector::k_end - Prospector::k_start);
    static const ui spacer_min =  21;
    static const ui spacer_max = 72;
    static const ui spacer_skip = (Prospector::spacer_min - 1);
    static const ui repeats_min = 3;
    static const ui size = 16;
    static const ui map_size_small = 64;
    static const ui map_size_big = 8000;
    static const ui bits = 2;
    static const ui repeat_tolerance_ratio = 8;

    struct Encoding {
        ui* encoding;
        ui* encoding_d;
        ui size;
    };

    void device_init();
    Prospector::Encoding get_genome_encoding(const char* genome, const ui genome_size);
    uc* get_qmap_small(const ui* genome_encoding, const ui genome_encoding_size);
    uc* get_qmap_big(const ui* genome_encoding, const ui genome_encoding_size, const ui* queries, const ui queries_size);
}

