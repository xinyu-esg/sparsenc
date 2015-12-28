#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "snc.h"

char usage[] = "usage: ./programName code_t dec_t sched_t datasize pcrate size_b size_g size_p [bufsize]\n\
                code_t  - code type: RAND, BAND, WINDWRAP\n\
                dec_t   - decoder type: GG, OA, BD, CBD\n\
                sched_t - scheduling type: TRIV, RAND, MLPI\n";
int main(int argc, char *argv[])
{
    if (argc != 9 && argc != 10) {
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

    int decoder_t;
    if (strcmp(argv[2], "GG") == 0)
        decoder_t = GG_DECODER;
    else if (strcmp(argv[2], "OA") == 0)
        decoder_t = OA_DECODER;
    else if (strcmp(argv[2], "BD") == 0)
        decoder_t = BD_DECODER;
    else if (strcmp(argv[2], "CBD") == 0)
        decoder_t = CBD_DECODER;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    int sched_t;
    if (strcmp(argv[3], "TRIV") == 0)
        sched_t = TRIV_SCHED;
    else if (strcmp(argv[3], "RAND") == 0)
        sched_t = RAND_SCHED;
    else if (strcmp(argv[3], "MLPI") == 0)
        sched_t = MLPI_SCHED;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    sp.datasize = atoi(argv[4]);
    sp.pcrate   = atof(argv[5]);
    sp.size_b   = atoi(argv[6]);
    sp.size_g   = atoi(argv[7]);
    sp.size_p   = atoi(argv[8]);
    int bufsize = 2;    // SNC buffer size, default is 2
    if (argc == 10)
        bufsize = atoi(argv[9]);
    sp.bpc      = 0;
    sp.bnc      = 0;

    srand( (int) time(0) );
    char *buf = malloc(sp.datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, sp.datasize);
    close(rnd);

    struct snc_context *sc;
    /* Create GNC encoding context */
    if ((sc = snc_create_enc_context(buf, sp)) == NULL) {
        fprintf(stderr, "Cannot create snc_context.\n");
        exit(1);
    }

    /* Create recoder buffer */
    struct snc_buffer *buffer;
    struct snc_metainfo *meta = snc_get_metainfo(sc);
    if ((buffer = snc_create_buffer(*meta, bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer.\n");
        exit(1);
    }

    /* Create decoder */
    struct snc_decoder *decoder = snc_create_decoder(sp, decoder_t);
    clock_t start, stop, dtime = 0;
    while (!snc_decoder_finished(decoder)) {
        struct snc_packet *pkt = snc_generate_packet(sc);
        snc_buffer_packet(buffer, pkt);

        struct snc_packet *rpkt = snc_recode_packet(buffer, sched_t);
        if (rpkt == NULL)
            continue;
        /* Measure decoding time */
        start = clock();
        snc_process_packet(decoder, rpkt);
        stop = clock();
        dtime += stop - start;
    }
    printf("dec-time: %.2f bufsize: %d ", ((double) dtime)/CLOCKS_PER_SEC, bufsize);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0) 
        fprintf(stderr, "recovered is NOT identical to original.\n");
    print_code_summary(dsc, snc_code_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    snc_free_buffer(buffer);
    snc_free_decoder(decoder);
    return 0;
}
