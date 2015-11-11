#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "slncEncoder.h"
#include "slncGGDecoder.h"

char usage[] = "usage: ./test.GGdecoder datasize size_b size_g size_p";
int main(int argc, char *argv[])
{
    if (argc != 6) {
        printf("%s\n", usage);
        exit(1);
    }
    long datasize = atoi(argv[1]);
    struct slnc_parameter sp;
    sp.pcrate   = atof(argv[2]);
    sp.size_b   = atoi(argv[3]);
    sp.size_g   = atoi(argv[4]);
    sp.size_p   = atoi(argv[5]);
    sp.type 	= BAND_GNC_CODE;

    srand( (int) time(0) );
    char *buf = malloc(datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, datasize);
    close(rnd);
    struct slnc_context *sc;

    if (create_slnc_context(buf, datasize, &sc, sp) != 0) {
        printf("Cannot create File Context.\n");
        return 1;
    }

    struct decoding_context_GG *dec_ctx = malloc(sizeof(struct decoding_context_GG));
    create_decoding_context_GG(dec_ctx, sc->meta.datasize, sp);
    while (dec_ctx->finished != 1) {
        struct coded_packet *pkt = generate_slnc_packet(sc);
        process_packet_GG(dec_ctx, pkt);
    }

    unsigned char *rec_buf = recover_data(dec_ctx->sc);
    if (memcmp(buf, rec_buf, datasize) != 0) 
        printf("recovered is NOT identical to original.\n");
    else
        printf("recovered is identical to original.\n");

    print_code_summary(&dec_ctx->sc->meta, dec_ctx->overhead, dec_ctx->operations);

    free_slnc_context(sc);
    free_decoding_context_GG(dec_ctx);
    return 0;
}
