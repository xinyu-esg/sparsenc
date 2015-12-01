#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "slncEncoder.h"
#include "slncDecoder.h"

char usage[] = "usage: ./programName code_t dec_t datasize pcrate size_b size_g size_p filename\n\
                       code_t - RAND, BAND, WINDWRAP\n\
                       dec_t  - GG, OA, BD, CBD\n";
int main(int argc, char *argv[])
{
    if (argc != 9) {
        printf("%s\n", usage);
        exit(1);
    }
    struct slnc_parameter sp;
    if (strcmp(argv[1], "RAND") == 0)
        sp.type = RAND_SLNC;
    else if (strcmp(argv[1], "BAND") == 0)
        sp.type = BAND_SLNC;
    else if (strcmp(argv[1], "WINDWRAP") == 0)
        sp.type = WINDWRAP_SLNC;
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
    long datasize = atoi(argv[3]);
    sp.pcrate   = atof(argv[4]);
    sp.size_b   = atoi(argv[5]);
    sp.size_g   = atoi(argv[6]);
    sp.size_p   = atoi(argv[7]);
    char *filename = argv[8];
    char *copyname = calloc(strlen(argv[8])+strlen(".dec.copy")+1, sizeof(char));
    strcat(copyname, filename);
    strcat(copyname, ".dec.copy");

    srand( (int) time(0) );
    /* Open file */
    FILE *fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("main: fopen %s failed\n", filename);
        exit(1);
    }
    fseek(fp, 0, SEEK_END);
    long filesize = ftell(fp);      // Determine filesize
    fseek(fp, 0, SEEK_SET);         // Seek back to start of the file

    int chunks = filesize % datasize == 0 ? filesize/datasize : filesize/datasize+1;
    printf("File splitted to %d chunks\n", chunks);
    struct slnc_context *sc;
    long remaining = filesize;
    long start = 0;
    while (chunks > 0) {
        /*
         * As a rough example, we create and free enc/dec context 
         * for each chunk. It's a waste of (alloc/init/free) time
         * to not reuse existing context. But we don't care because
         * this is just a test of APIs.
         */
        long toEncode = remaining > datasize ? datasize : remaining; 
        if (slnc_create_enc_context(NULL, toEncode, &sc, sp) != 0) {
            fprintf(stderr, "Cannot create File Context.\n");
            return 1;
        }
        slnc_load_file_to_context(fp, start, sc);
        remaining -= toEncode;
        start += toEncode;
        chunks--;

        slnc_create_dec_context(decoder_type, sc->meta.datasize, sp);
        clock_t start, stop, dtime = 0;
        while (slnc_is_dec_finished(decoder_type) != 1) {
            struct slnc_packet *pkt = slnc_generate_packet(sc);
            /* Measure decoding time */
            start = clock();
            slnc_process_packet(decoder_type, pkt);
            stop = clock();
            dtime += stop - start;
        }
        printf("chunk %d || dec-time: %.2f ", chunks, ((double) dtime)/CLOCKS_PER_SEC);

        struct slnc_context *dsc = slnc_decoded_context(decoder_type);
        FILE *wfp = fopen(copyname, "a");
        if (wfp == NULL) {
            printf("main: fopen %s failed.\n", copyname);
            exit(1);
        }
        slnc_recover_to_file(wfp, dsc);
        fclose(wfp);

        print_code_summary(&(dsc->meta), slnc_dec_overhead(decoder_type), slnc_dec_operations(decoder_type));
        slnc_free_enc_context(sc);
        slnc_free_dec_context(decoder_type);
    }


    return 0;
}
