#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "snc.h"

char usage[] = "usage: ./programName code_t dec_t datasize pcrate size_b size_g size_p\n\
                       code_t - RAND, BAND, WINDWRAP\n\
                       dec_t  - GG, OA, BD, CBD\n";
int main(int argc, char *argv[])
{
    if (argc != 8) {
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

    srand( (int) time(0) );
    char *buf = malloc(sp.datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, sp.datasize);
    close(rnd);

    struct snc_context *sc;
    if ((sc = snc_create_enc_context(buf, sp)) == NULL) {
        fprintf(stderr, "Cannot create File Context.\n");
        return 1;
    }

    struct snc_decoder *decoder = snc_create_decoder(sp, decoder_type);
    clock_t start, stop, dtime = 0;
    while (snc_decoder_finished(decoder) != 1) {
        struct snc_packet *pkt = snc_generate_packet(sc);
        /* Measure decoding time */
        start = clock();
		snc_process_packet(decoder, pkt);
        stop = clock();
        dtime += stop - start;
    }
    printf("dec-time: %.2f ", ((double) dtime)/CLOCKS_PER_SEC);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0) 
        fprintf(stderr, "recovered is NOT identical to original.\n");

    print_code_summary(dsc, snc_code_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    snc_free_decoder(decoder);
    return 0;
}
