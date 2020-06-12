#pragma once

typedef unsigned long long ull;
typedef unsigned int ui;
typedef signed long int li;
typedef unsigned char uc;

namespace Prospector
{
    static const ui k_start = 20;
    static const ui k_end = 60;
    static const ui k_count = (k_end - k_start);
    static const ui spacer_min =  21;
    static const ui spacer_max = 72;
    static const ui spacer_skip = (spacer_min - 1);
    static const ui repeats_min = 3;
    static const ui size = 16;
    static const ui map_size_small = 64;
    static const ui bits = 2;
    static const ui repeat_tolerance_ratio = 15;

    struct Encoding {
        ui* h;
        ui* d;
        ui size;
        ui bytes;
    };

    void device_init();
    Encoding get_genome_encoding(const char*, const ui);
    uc* get_qmap_small(const ui*, const ui);
}


