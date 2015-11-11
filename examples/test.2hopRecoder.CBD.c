#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "slncEncoder.h"
#include "slncRecoder.h"
#include "slncCBDDecoder.h"

char usage[] = "usage: ./programName datasize pcrate size_b size_g size_p [bufsize]";
int main(int argc, char *argv[])
{
    if (argc != 6 && argc != 7) {
        printf("%s\n", usage);
        exit(1);
    }
    long datasize = atoi(argv[1]);
    struct slnc_parameter sp;
    sp.pcrate   = atof(argv[2]);
    sp.size_b   = atoi(argv[3]);
    sp.size_g   = atoi(argv[4]);
    sp.size_p   = atoi(argv[5]);
    sp.type     = BAND_GNC_CODE;
    int bufsize = 2;    // Buffer size of each generation, default is 2
    if (argc == 7)
        bufsize = atoi(argv[6]);


    srand( (int) time(0) );
    char *buf = malloc(datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, datasize);
    close(rnd);
    struct slnc_context *sc;

    /* Create GNC encoding context */
    if (create_slnc_context(buf, datasize, &sc, sp) != 0) {
        fprintf(stderr, "Cannot create File Context.\n");
        return 1;
    }

    /* Create recoding context */
    struct slnc_recoding_context *rc = malloc(sizeof(struct slnc_recoding_context));
    if (create_recoding_context(rc, sc->meta, bufsize) != 0) {
        fprintf(stderr, "Cannot create recoding context.\n");
        return 1;
    }

    /* Create CBD decoding context */
    struct decoding_context_CBD *dec_ctx = malloc(sizeof(struct decoding_context_CBD));
    create_decoding_context_CBD(dec_ctx, sc->meta.datasize, sp);
    clock_t start, stop, dtime = 0;
    while (dec_ctx->finished != 1) {
        struct coded_packet *pkt = generate_slnc_packet(sc);
        buffer_packet(rc, pkt);

        struct coded_packet *rpkt = generate_recoded_packet(rc, MLPI_SCHED);
        if (rpkt == NULL)
            continue;
        /* Measure decoding time */
        start = clock();
        process_packet_CBD(dec_ctx, rpkt);
        stop = clock();
        dtime += stop - start;
    }
    printf("dec-time: %.2f bufsize: %d ", ((double) dtime)/CLOCKS_PER_SEC, bufsize);

    unsigned char *rec_buf = recover_data(dec_ctx->sc);
    if (memcmp(buf, rec_buf, datasize) != 0) 
        fprintf(stderr, "recovered is NOT identical to original.\n");

    print_code_summary(&dec_ctx->sc->meta, dec_ctx->overhead, dec_ctx->operations);

    free_slnc_context(sc);
    free_recoding_buffer(rc);
    free(rc);
    free_decoding_context_CBD(dec_ctx);
    return 0;
}
