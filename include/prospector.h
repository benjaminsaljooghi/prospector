#pragma once

typedef long long ll;
typedef unsigned int ui;
typedef unsigned int kmer;
typedef unsigned long long ull;
typedef unsigned char uc;

namespace Prospector
{
    static const ui k_start = 20;
    static const ui k_end = 60;
    static const ull k_count = (k_end - k_start);
    static const ull spacer_min =  21;
    static const ull spacer_max = 72;
    static const ull spacer_skip = (spacer_min - 1);
    static const ull repeats_min = 3;
    static const ull size = 16;
    static const ull map_size_small = 64;
    static const ull bits = 2;
    static const ull repeat_tolerance_ratio = 10; // lower is more sensitive
    static const ull repeat_tolerance_ratio_sensitive = 4; // lower is more sensitive

    struct Encoding {
        kmer* h;
        kmer* d;
        ull size;
        ull bytes;
    };

    void device_init();
    Encoding get_genome_encoding(const char*, const ull);
    uc* get_qmap_small(const kmer* encoding, const ull encoding_size);

    void free_encoding(Encoding);
}


