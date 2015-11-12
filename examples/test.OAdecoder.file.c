#include <unistd.h>
#include <fcntl.h>
#include "common.h"
#include "slncEncoder.h"
#include "slncOADecoder.h"

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
    int slnc_type   = BAND_SLNC;

    //char filename[] = "test.file";
    srand( (int) time(0) );
    struct slnc_context *sc;
    /* Open file */
    FILE *fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("main: fopen %s failed\n", filename);
        exit(1);
    }

    if (slnc_create_enc_context_from_file(fp, &sc, size_b, size_g, size_p, slnc_type) != 0) {
        printf("Cannot create File Context.\n");
        return 1;
    }
    fclose(fp);

    printf("File size: %ld\n", sc->meta.datasize);
    printf("Number of packets: %d\n", sc->meta.snum);

    struct slnc_dec_context_OA *dec_ctx = malloc(sizeof(struct slnc_dec_context_OA));
    slnc_create_dec_context_OA(dec_ctx, sc->meta.datasize, sc->meta.size_b, sc->meta.size_g, sc->meta.size_p, sc->meta.type, 0);
    while (!dec_ctx->finished) {
        struct slnc_packet *pkt = slnc_generate_packet(sc);
        slnc_process_packet_OA(dec_ctx, pkt);
    }
    printf("overhead: %f computation: %f/symbol\n", 
            (double) dec_ctx->overhead/dec_ctx->sc->meta.snum,
            (double) dec_ctx->operations/dec_ctx->sc->meta.snum/(dec_ctx->sc->meta.size_p));

    char *copyname = calloc(strlen(argv[1])+strlen(".dec.copy")+1, sizeof(char));
    strcat(copyname, filename);
    strcat(copyname, ".dec.copy");
    FILE *wfp = fopen(copyname, "a");
    if (wfp == NULL) {
        printf("main: fopen %s failed.\n", copyname);
        exit(1);
    }
    slnc_recover_data_to_file(wfp, dec_ctx->sc);
    fclose(wfp);
    slnc_free_enc_context(sc);
    slnc_free_dec_context_OA(dec_ctx);
    return 0;
}
