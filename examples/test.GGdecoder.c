#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gncEncoder.h"
#include "gncGGDecoder.h"

char usage[] = "usage: ./test.GGdecoder datasize size_b size_g size_p";
int main(int argc, char *argv[])
{
    if (argc != 6) {
        printf("%s\n", usage);
        exit(1);
    }
    long datasize = atoi(argv[1]);
    struct gnc_parameter gp;
    gp.pcrate   = atof(argv[2]);
    gp.size_b   = atoi(argv[3]);
    gp.size_g   = atoi(argv[4]);
    gp.size_p   = atoi(argv[5]);
    gp.type 	= RAND_GNC_CODE;

    srand( (int) time(0) );
    char *buf = malloc(datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, datasize);
    close(rnd);
    struct gnc_context *gc;

    if (create_gnc_context(buf, datasize, &gc, gp) != 0) {
        printf("Cannot create File Context.\n");
        return 1;
    }

    struct decoding_context_GG *dec_ctx = malloc(sizeof(struct decoding_context_GG));
    create_decoding_context_GG(dec_ctx, gc->meta.datasize, gp);
    while (dec_ctx->finished != 1) {
        struct coded_packet *pkt = generate_gnc_packet(gc);
        process_packet_GG(dec_ctx, pkt);
    }

    unsigned char *rec_buf = recover_data(dec_ctx->gc);
    if (memcmp(buf, rec_buf, datasize) != 0) 
        printf("recovered is NOT identical to original.\n");
    else
        printf("recovered is identical to original.\n");

    print_code_summary(&dec_ctx->gc->meta, dec_ctx->overhead, dec_ctx->operations);

    free_gnc_context(gc);
    free_decoding_context_GG(dec_ctx);
    return 0;
}
