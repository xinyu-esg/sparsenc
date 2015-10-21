#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gncEncoder.h"
#include "gncBandDecoder.h"

char usage[] = "usage: ./programName datasize pcrate size_b size_g size_p";
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
	gp.type     = BAND_GNC_CODE;

	srand( (int) time(0) );
	char *buf = malloc(datasize);
	int rnd=open("/dev/urandom", O_RDONLY);
	read(rnd, buf, datasize);
	close(rnd);
	struct gnc_context *gc;

	// Construct a GNC code (32, 40, 1024)
	if (create_gnc_context(buf, datasize, &gc, gp) != 0) {
		fprintf(stderr, "Cannot create File Context.\n");
		return 1;
	}

	struct decoding_context_BD *dec_ctx = malloc(sizeof(struct decoding_context_BD));
	create_decoding_context_BD(dec_ctx, gc->meta.datasize, gp);
	clock_t start, stop, dtime = 0;
	while (dec_ctx->finished != 1) {
		struct coded_packet *pkt = generate_gnc_packet(gc);
		/* Measure decoding time */
		start = clock();
		process_packet_BD(dec_ctx, pkt);
		stop = clock();
		dtime += stop - start;
	}
	printf("dec-time: %.2f ", ((double) dtime)/CLOCKS_PER_SEC);

	unsigned char *rec_buf = recover_data(dec_ctx->gc);
	if (memcmp(buf, rec_buf, datasize) != 0) 
		fprintf(stderr, "recovered is NOT identical to original.\n");

	print_code_summary(&dec_ctx->gc->meta, dec_ctx->overhead, dec_ctx->operations);

	free_gnc_context(gc);
	free_decoding_context_BD(dec_ctx);
	return 0;
}
