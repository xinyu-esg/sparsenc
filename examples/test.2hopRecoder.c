#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "usage: ./programName code_t dec_t sched_t datasize pcrate size_b size_g size_p bpc bnc sys pe1 pe2 [bufsize]\n\
                code_t  - code type: RAND, BAND, WINDWRAP\n\
                dec_t   - decoder type: GG, OA, BD, CBD\n\
                sched_t - scheduling type: TRIV, RAND, MLPI, NURAND\n\
                bpc      - Use binary precode (0 or 1)\n\
                bnc      - Use binary network code (0 or 1)\n\
                sys      - Systematic code (0 or 1)\n\
                pe1, pe2 - erasure rates of the two hops\n";
int main(int argc, char *argv[])
{
    if (argc != 14 && argc != 15) {
        printf("%s\n", usage);
        exit(1);
    }
    struct snc_parameters sp;
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
    else if (strcmp(argv[3], "NURAND") == 0)
        sched_t = NURAND_SCHED;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    sp.datasize = atoi(argv[4]);
    sp.pcrate   = atof(argv[5]);
    sp.size_b   = atoi(argv[6]);
    sp.size_g   = atoi(argv[7]);
    sp.size_p   = atoi(argv[8]);
    sp.bpc      = atoi(argv[9]);
    sp.bnc      = atoi(argv[10]);
    sp.sys      = atoi(argv[11]);
    sp.seed     = -1;
    double pe1  = atof(argv[12]);
    double pe2  = atof(argv[13]);
    int bufsize = 2;    // SNC buffer size, default is 2
    if (argc == 15)
        bufsize = atoi(argv[14]);

    char *ur = getenv("SNC_NONUNIFORM_RAND");
    if ( ur != NULL && atoi(ur) == 1) {
        if (sp.type != BAND_SNC || sp.size_b != 1) {
            printf("Non-Uniform Random Scheduling can only be used for BAND code with size_b=1.\n");
            exit(1);
        }
    }

    if (sched_t == NURAND_SCHED) {
        if (sp.type != BAND_SNC || sp.size_b != 1) {
            printf("Non-Uniform Random Scheduling can only be used for BAND code with size_b=1.\n");
            exit(1);
        }
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec
    char *buf = malloc(sp.datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, sp.datasize);
    close(rnd);

    struct snc_context *sc;
    /* Create GNC encoding context */
    if ((sc = snc_create_enc_context(buf, &sp)) == NULL) {
        fprintf(stderr, "Cannot create snc_context.\n");
        exit(1);
    }

    /* Create recoder buffer */
    struct snc_buffer *buffer;
    if ((buffer = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer.\n");
        exit(1);
    }

    /* Create decoder */
    sp.seed = (snc_get_parameters(sc))->seed;
    struct snc_decoder *decoder = snc_create_decoder(&sp, decoder_t);
    struct snc_packet *rpkt = snc_alloc_empty_packet(&sp);
    clock_t start, stop, dtime = 0;
    // double pe1 = 0.7;
    // double pe2 = 0.4;
    while (!snc_decoder_finished(decoder)) {
        struct snc_packet *pkt = snc_generate_packet(sc);
        if (rand() % 100 >= pe1 *100)
            snc_buffer_packet(buffer, pkt);
        else
            snc_free_packet(pkt);
        
        //struct snc_packet *rpkt = snc_recode_packet(buffer, sched_t);
        int ret = snc_recode_packet_im(buffer, rpkt, sched_t);
        if (ret == -1)
            continue;
        if (rand() % 100 < pe2 * 100)
            continue;
        /* Measure decoding time */
        start = clock();
        snc_process_packet(decoder, rpkt);
        stop = clock();
        dtime += stop - start;
    }
    printf("dec-time: %.2f bufsize: %d pe1: %.3f pe2: %.3f ", ((double) dtime)/CLOCKS_PER_SEC, bufsize, pe1, pe2);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0)
        fprintf(stderr, "recovered is NOT identical to original.\n");
    print_code_summary(dsc, snc_decode_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    snc_free_buffer(buffer);
    snc_free_packet(rpkt);
    snc_free_decoder(decoder);
    return 0;
}
