#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gncEncoder.h"
#include "gncRecoder.h"
#include "gncCBDDecoder.h"

char usage[] = "usage: ./programName datasize pcrate size_b size_g size_p [bufsize]";
int main(int argc, char *argv[])
{
	if (argc != 6 && argc != 7) {
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
	int bufsize = 2;    // Buffer size of each generation, default is 100
	if (argc == 7)
		bufsize = atoi(argv[6]);


	srand( (int) time(0) );
	char *buf = malloc(datasize);
	int rnd=open("/dev/urandom", O_RDONLY);
	read(rnd, buf, datasize);
	close(rnd);
	struct gnc_context *gc;

	/* Create GNC encoding context */
	if (create_gnc_context(buf, datasize, &gc, gp) != 0) {
		fprintf(stderr, "Cannot create File Context.\n");
		return 1;
	}

	/* Create recoding context */
	struct gnc_recoding_context *rc = malloc(sizeof(struct gnc_recoding_context));
    if (create_recoding_context(rc, gc->meta, bufsize) != 0) {
		fprintf(stderr, "Cannot create recoding context.\n");
		return 1;
	}

	/* Create CBD decoding context */
	struct decoding_context_CBD *dec_ctx = malloc(sizeof(struct decoding_context_CBD));
	create_decoding_context_CBD(dec_ctx, gc->meta.datasize, gp);
	clock_t start, stop, dtime = 0;
	while (dec_ctx->finished != 1) {
		struct coded_packet *pkt = generate_gnc_packet(gc);
		buffer_packet(rc, pkt);
		
		struct coded_packet *rpkt = generate_recoded_packet(rc, RAND_SCHED);
		if (rpkt == NULL)
			continue;
		/* Measure decoding time */
		start = clock();
		process_packet_CBD(dec_ctx, rpkt);
		stop = clock();
		dtime += stop - start;
	}
	printf("dec-time: %.2f ", ((double) dtime)/CLOCKS_PER_SEC);

	unsigned char *rec_buf = recover_data(dec_ctx->gc);
	if (memcmp(buf, rec_buf, datasize) != 0) 
		fprintf(stderr, "recovered is NOT identical to original.\n");

	print_code_summary(&dec_ctx->gc->meta, dec_ctx->overhead, dec_ctx->operations);

	free_gnc_context(gc);
	free_recoding_buffer(rc);
	free(rc);
	free_decoding_context_CBD(dec_ctx);
	return 0;
}
