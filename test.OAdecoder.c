#include <unistd.h>
#include <fcntl.h>
#include "common.h"
#include "gncEncoder.h"
#include "gncOADecoder.h"

int main(int argc, char *argv[])
{
	//char filename[] = "test.file";
	srand( (int) time(0) );
	long datasize = 1599800;
	unsigned char *buf = malloc(datasize);
	int rnd=open("/dev/urandom", O_RDONLY);
	read(rnd, buf, datasize);
	close(rnd);
	struct gnc_context *gc;

	// Construct a GNC code (32, 40, 1024)
	if (create_gnc_context(buf, datasize, &gc, 32, 42, 1024) != 0) {
		printf("Cannot create File Context.\n");
		return 1;
	}

	printf("File size: %ld\n", gc->meta.datasize);
	printf("Number of packets: %d\n", gc->meta.snum);

	struct decoding_context_OA *dec_ctx = malloc(sizeof(struct decoding_context_OA));
	create_decoding_context_OA(dec_ctx, datasize, 32, 42, 1024, 0);
	while (dec_ctx->finished != 1) {
		struct coded_packet *pkt = generate_gnc_packet(gc);
		process_packet_OA(dec_ctx, pkt);
	}
	printf("overhead: %f computation: %f/symbol\n", 
					(double) dec_ctx->overhead/dec_ctx->gc->meta.snum,
					(double) dec_ctx->operations/dec_ctx->gc->meta.snum/(dec_ctx->gc->meta.size_p));

	unsigned char *rec_buf = recover_data(dec_ctx->gc);
	if (memcmp(buf, rec_buf, datasize) != 0) 
		printf("recovered is NOT identical to original.\n");
	else
		printf("recovered is identical to original.\n");

	free_gnc_context(gc);
	//free_decoding_context_GG(dec_ctx);
	return 0;
}
