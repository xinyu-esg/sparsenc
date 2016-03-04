#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "usage: ./programName code_t dec_t datasize pcrate size_b size_g size_p filename\n\
                       code_t - RAND, BAND, WINDWRAP\n\
                       dec_t  - GG, OA, BD, CBD\n";
int main(int argc, char *argv[])
{
    if (argc != 9) {
        printf("%s\n", usage);
        exit(1);
    }
    struct snc_parameter sp;
    if (strcmp(argv[1], "RAND") == 0)
        sp.type = RAND_SNC;
    else if (strcmp(argv[1], "BAND") == 0)
        sp.type = BAND_SNC;
    else if (strcmp(argv[1], "WINDWRAP") == 0)
        sp.type = WINDWRAP_SNC;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    int decoder_type;
    if (strcmp(argv[2], "GG") == 0)
        decoder_type = GG_DECODER;
    else if (strcmp(argv[2], "OA") == 0)
        decoder_type = OA_DECODER;
    else if (strcmp(argv[2], "BD") == 0)
        decoder_type = BD_DECODER;
    else if (strcmp(argv[2], "CBD") == 0)
        decoder_type = CBD_DECODER;
    else {
        printf("%s\n", usage);
        exit(1);
    }
    sp.datasize = atoi(argv[3]);
    sp.pcrate   = atof(argv[4]);
    sp.size_b   = atoi(argv[5]);
    sp.size_g   = atoi(argv[6]);
    sp.size_p   = atoi(argv[7]);
    sp.bpc      = 0;
    sp.bnc      = 0;
    sp.sys      = atoi(argv[10]);
    sp.seed     = -1;  // Initialize seed as -1
    char *filename = argv[8];
    char *copyname = calloc(strlen(argv[8])+strlen(".dec.copy")+1, sizeof(char));
    strcat(copyname, filename);
    strcat(copyname, ".dec.copy");

    srand( (int) time(0) );
    // Get filesize in bytes
    struct stat sb;
    if (stat(filename, &sb) < -1) {
        printf("main: cannot obtain %s's stat\n", filename);
        exit(1);
    }
    long filesize = sb.st_size;

    int chunks = filesize % sp.datasize == 0 ? filesize/sp.datasize : filesize/sp.datasize+1;
    printf("File size: %d splitted to %d chunks\n", filesize, chunks);
    struct snc_context *sc;
    long remaining = filesize;
    long start = 0;
    while (chunks > 0) {
        /*
         * As a rough example, we create and free enc/dec context 
         * for each chunk. It's a waste of (alloc/init/free) time
         * to not reuse existing context. But we don't care because
         * this is just a test of APIs.
         */
        long toEncode = remaining > sp.datasize ? sp.datasize : remaining; 
        sp.datasize = toEncode;
        if ((sc = snc_create_enc_context(NULL, &sp)) == NULL) {
            fprintf(stderr, "Cannot create File Context.\n");
            return 1;
        }
        snc_load_file_to_context(filename, start, sc);
        remaining -= toEncode;
        start += toEncode;
        chunks--;

        struct snc_decoder *decoder = snc_create_decoder(&sp, decoder_type);
        clock_t start, stop, dtime = 0;
        while (snc_decoder_finished(decoder) != 1) {
            struct snc_packet *pkt = snc_generate_packet(sc);
            /* Measure decoding time */
            start = clock();
            snc_process_packet(decoder, pkt);
            stop = clock();
            dtime += stop - start;
        }
        printf("chunk %d || dec-time: %.2f ", chunks, ((double) dtime)/CLOCKS_PER_SEC);

        struct snc_context *dsc = snc_get_enc_context(decoder);
        snc_recover_to_file(copyname, dsc);

        print_code_summary(dsc, snc_code_overhead(decoder), snc_decode_cost(decoder));

        snc_free_enc_context(sc);
        snc_free_decoder(decoder);
    }


    return 0;
}
