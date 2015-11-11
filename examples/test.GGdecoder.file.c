#include <unistd.h>
#include <fcntl.h>
#include "slncEncoder.h"
#include "slncGGDecoder.h"

/*
 * This binary tests GG (generation-by-generation) decoder on GNC code
 * encoded from a file.
 *
 */

char usage[] = "usage: ./test.GGdecoder.file filename size_b size_g size_p";

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
    int slnc_type   = RAND_GNC_CODE;

    //char filename[] = "test.file";
    srand( (int) time(0) );
    struct slnc_context *sc;
    /* Open file */
    FILE *fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("main: fopen %s failed\n", filename);
        exit(1);
    }

    if (create_slnc_context_from_file(fp, &sc, size_b, size_g, size_p, slnc_type) != 0) {
        printf("Cannot create File Context.\n");
        return 1;
    }
    if (load_file_to_slnc_context(fp, sc) != 0) {
        printf("Load file to slnc_context failed.\n");
        return 1;
    }
    fclose(fp);

    printf("File size: %ld\n", sc->meta.datasize);
    printf("Number of packets: %d\n", sc->meta.snum);

    struct decoding_context_GG *dec_ctx = malloc(sizeof(struct decoding_context_GG));
    create_decoding_context_GG(dec_ctx, sc->meta.datasize, sc->meta.size_b, sc->meta.size_g, sc->meta.size_p, sc->meta.type);
    while (!dec_ctx->finished) {
        struct coded_packet *pkt = generate_slnc_packet(sc);
        process_packet_GG(dec_ctx, pkt);
    }

    char *copyname = calloc(strlen(argv[1])+strlen(".dec.copy")+1, sizeof(char));
    strcat(copyname, filename);
    strcat(copyname, ".dec.copy");
    FILE *wfp = fopen(copyname, "a");
    if (wfp == NULL) {
        printf("main: fopen %s failed.\n", copyname);
        exit(1);
    }
    recover_data_to_file(wfp, dec_ctx->sc);
    fclose(wfp);
    print_code_summary(&dec_ctx->sc->meta, dec_ctx->overhead, dec_ctx->operations);
    free_slnc_context(sc);
    free_decoding_context_GG(dec_ctx);
    return 0;
}
