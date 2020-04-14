// #include "stdafx.h"
// #include "crispr.h"
// #include "util.h"

#define SIZE 16
#define MAP_SIZE 64
#define BITS 2
#define MUTANT_TOLERANCE_RATIO 6


#define K_START 20
#define K_END 60
#define K_COUNT (K_END - K_START)


#define REPEAT_MIN 20
#define REPEAT_MAX 60
#define SPACER_MIN 21
#define SPACER_MAX 72
#define SPACER_SKIP (SPACER_MIN - 1)
#define REPEATS_MIN 3
#define SCAN_DOMAIN SPACER_MAX
#define ALLOW_DISCREPANT_LENGTHS false
#define MIN_REPEATS 3
//#define CRISPR_BUFFER 50
#define printf_BYTE_FORMAT_ALIGN 10


typedef unsigned long long ull;
typedef unsigned int ui;
typedef signed long int li;
typedef unsigned char uc;


namespace Prospector
{
    struct Encoding {
        ui* encoding;
        ui* d_encoding;
        ui size;
    };

    void device_init();
    Encoding get_genome_encoding(const char* genome, ui genome_size);

    uc* get_qmap_small(ui* genome_encoding, ui genome_encoding_size);
    uc* get_qmap_big(const ui* genome_encoding, const ui genome_encoding_size, const ui* queries, const ui queries_size, ui map_size);
}

