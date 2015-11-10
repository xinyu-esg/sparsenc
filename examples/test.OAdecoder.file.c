#include <unistd.h>
#include <fcntl.h>
#include "common.h"
#include "gncEncoder.h"
#include "gncOADecoder.h"

/*
 * This binary tests GG (generation-by-generation) decoder on GNC code
 * encoded from a file.
 *
 */

char usage[] = "usage: ./test.OAdecoder.file filename size_b size_g size_p";

int main(int argc, char *argv[])
{
    if (argc != 5) {
        printf("%s\n", usage);
        exit(1);
    }

    char *filename = argv[1];
    int size_b     = atoi(argv[2]);
    int size_g     = atoi(argv[3]);
    int size_p     = atoi(argv[4]);
    int gnc_type   = BAND_GNC_CODE;

    //char filename[] = "test.file";
    srand( (int) time(0) );
    struct gnc_context *gc;
    /* Open file */
    FILE *fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("main: fopen %s failed\n", filename);
        exit(1);
    }

    if (create_gnc_context_from_file(fp, &gc, size_b, size_g, size_p, gnc_type) != 0) {
        printf("Cannot create File Context.\n");
        return 1;
    }
    fclose(fp);

    printf("File size: %ld\n", gc->meta.datasize);
    printf("Number of packets: %d\n", gc->meta.snum);

    struct decoding_context_OA *dec_ctx = malloc(sizeof(struct decoding_context_OA));
    create_decoding_context_OA(dec_ctx, gc->meta.datasize, gc->meta.size_b, gc->meta.size_g, gc->meta.size_p, gc->meta.type, 0);
    while (!dec_ctx->finished) {
        struct coded_packet *pkt = generate_gnc_packet(gc);
        process_packet_OA(dec_ctx, pkt);
    }
    printf("overhead: %f computation: %f/symbol\n", 
            (double) dec_ctx->overhead/dec_ctx->gc->meta.snum,
            (double) dec_ctx->operations/dec_ctx->gc->meta.snum/(dec_ctx->gc->meta.size_p));

    char *copyname = calloc(strlen(argv[1])+strlen(".dec.copy")+1, sizeof(char));
    strcat(copyname, filename);
    strcat(copyname, ".dec.copy");
    FILE *wfp = fopen(copyname, "a");
    if (wfp == NULL) {
        printf("main: fopen %s failed.\n", copyname);
        exit(1);
    }
    recover_data_to_file(wfp, dec_ctx->gc);
    fclose(wfp);
    free_gnc_context(gc);
    free_decoding_context_OA(dec_ctx);
    return 0;
}
