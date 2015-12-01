#include <stdlib.h>
#include "decoderGG.h"
#include "decoderOA.h"
#include "decoderBD.h"
#include "decoderCBD.h"
#include "slncDecoder.h"

/**
 * decoding context as a static global variable
 **/
static void *dec_ctx = NULL;

int slnc_create_dec_context(int d_type, long datasize, struct slnc_parameter sp) {
	if (dec_ctx != NULL) {
		printf("Fatal: decoding context is allocated. Something must be wrong.\n");
		exit(1);
	}
	switch (d_type) {
	case GG_DECODER:
        dec_ctx = malloc(sizeof(struct decoding_context_GG));
        create_dec_context_GG(((struct decoding_context_GG *) dec_ctx), datasize, sp);
		break;
	case OA_DECODER:
        dec_ctx = malloc(sizeof(struct decoding_context_OA));
        create_dec_context_OA(((struct decoding_context_OA *) dec_ctx), datasize, sp, 0);
		break;
	case BD_DECODER:
        dec_ctx = malloc(sizeof(struct decoding_context_BD));
        create_dec_context_BD(((struct decoding_context_BD *) dec_ctx), datasize, sp);
		break;
	case CBD_DECODER:
        dec_ctx = malloc(sizeof(struct decoding_context_CBD));
        create_dec_context_CBD(((struct decoding_context_CBD *) dec_ctx), datasize, sp);
		break;
	}
	return (0);
}

void slnc_process_packet(int d_type, struct slnc_packet *pkt) {
	switch (d_type) {
	case GG_DECODER:
		process_packet_GG(((struct decoding_context_GG *) dec_ctx), pkt);
		break;
	case OA_DECODER:
		process_packet_OA(((struct decoding_context_OA *) dec_ctx), pkt);
		break;
	case BD_DECODER:
		process_packet_BD(((struct decoding_context_BD *) dec_ctx), pkt);
		break;
	case CBD_DECODER:
		process_packet_CBD(((struct decoding_context_CBD *) dec_ctx), pkt);
		break;
	}
	return;
}

int slnc_is_dec_finished(int d_type) {
	switch (d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) dec_ctx)->finished;
	case OA_DECODER:
		return ((struct decoding_context_OA *) dec_ctx)->finished;
	case BD_DECODER:
		return ((struct decoding_context_BD *) dec_ctx)->finished;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) dec_ctx)->finished;
	}
	return 0;
}

struct slnc_context *slnc_decoded_context(int d_type) {
	switch (d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) dec_ctx)->sc;
	case OA_DECODER:
		return ((struct decoding_context_OA *) dec_ctx)->sc;
	case BD_DECODER:
		return ((struct decoding_context_BD *) dec_ctx)->sc;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) dec_ctx)->sc;
	}
	return 0;
}

int slnc_dec_overhead(int d_type) {
	switch (d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) dec_ctx)->overhead;
	case OA_DECODER:
		return ((struct decoding_context_OA *) dec_ctx)->overhead;
	case BD_DECODER:
		return ((struct decoding_context_BD *) dec_ctx)->overhead;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) dec_ctx)->overhead;
	}
	return 0;
}

long long slnc_dec_operations(int d_type) {
	switch (d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) dec_ctx)->operations;
	case OA_DECODER:
		return ((struct decoding_context_OA *) dec_ctx)->operations;
	case BD_DECODER:
		return ((struct decoding_context_BD *) dec_ctx)->operations;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) dec_ctx)->operations;
	}
	return 0;
}

void slnc_free_dec_context(int d_type) {
	switch (d_type) {
	case GG_DECODER:
        free_dec_context_GG(((struct decoding_context_GG *) dec_ctx));
		break;
	case OA_DECODER:
        free_dec_context_OA(((struct decoding_context_OA *) dec_ctx));
		break;
	case BD_DECODER:
        free_dec_context_BD(((struct decoding_context_BD *) dec_ctx));
		break;
	case CBD_DECODER:
        free_dec_context_CBD(((struct decoding_context_CBD *) dec_ctx));
		break;
	}
	dec_ctx = NULL;
	return;
}
