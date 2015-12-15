#include <stdlib.h>
#include "decoderGG.h"
#include "decoderOA.h"
#include "decoderBD.h"
#include "decoderCBD.h"
#include "snc.h"

/* Definition of SNC decoder
 * It includes an allocated decoder context which contains encode
 * information (i.e., encode context) and miscellaneous structs that
 * are used during decoding. Different decoder types have different 
 * decode context, so d_type is used as an indicator.
 */
struct snc_decoder {
    void   *dec_ctx;        // decoder context
    int    d_type;          // decoder type
};

struct snc_decoder *snc_create_decoder(struct snc_parameter sp, int d_type) {
    struct snc_decoder *decoder = malloc(sizeof(struct snc_decoder));
    if (decoder == NULL)
        return NULL;

    decoder->d_type = d_type;

	switch (d_type) {
	case GG_DECODER:
        decoder->dec_ctx = malloc(sizeof(struct decoding_context_GG));
        if (decoder->dec_ctx == NULL) 
            goto failure; 
        create_dec_context_GG(((struct decoding_context_GG *) decoder->dec_ctx), sp);
		break;
	case OA_DECODER:
        decoder->dec_ctx = malloc(sizeof(struct decoding_context_OA));
        if (decoder->dec_ctx == NULL) 
            goto failure; 
        create_dec_context_OA(((struct decoding_context_OA *) decoder->dec_ctx), sp, 0);
		break;
	case BD_DECODER:
        decoder->dec_ctx = malloc(sizeof(struct decoding_context_BD));
        if (decoder->dec_ctx == NULL) 
            goto failure; 
        create_dec_context_BD(((struct decoding_context_BD *) decoder->dec_ctx), sp);
		break;
	case CBD_DECODER:
        decoder->dec_ctx = malloc(sizeof(struct decoding_context_CBD));
        if (decoder->dec_ctx == NULL) 
            goto failure; 
        create_dec_context_CBD(((struct decoding_context_CBD *) decoder->dec_ctx), sp);
		break;
	}
	return decoder; 
failure:
    free(decoder);
    return NULL;
}

void snc_process_packet(struct snc_decoder *decoder,  struct snc_packet *pkt) {
	switch (decoder->d_type) {
	case GG_DECODER:
		process_packet_GG(((struct decoding_context_GG *) decoder->dec_ctx), pkt);
		break;
	case OA_DECODER:
		process_packet_OA(((struct decoding_context_OA *) decoder->dec_ctx), pkt);
		break;
	case BD_DECODER:
		process_packet_BD(((struct decoding_context_BD *) decoder->dec_ctx), pkt);
		break;
	case CBD_DECODER:
		process_packet_CBD(((struct decoding_context_CBD *) decoder->dec_ctx), pkt);
		break;
	}
	return;
}

int snc_decoder_finished(struct snc_decoder *decoder) {
	switch (decoder->d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) decoder->dec_ctx)->finished;
	case OA_DECODER:
		return ((struct decoding_context_OA *) decoder->dec_ctx)->finished;
	case BD_DECODER:
		return ((struct decoding_context_BD *) decoder->dec_ctx)->finished;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) decoder->dec_ctx)->finished;
	}
	return 0;
}

struct snc_context *snc_get_enc_context(struct snc_decoder *decoder) {
	switch (decoder->d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) decoder->dec_ctx)->sc;
	case OA_DECODER:
		return ((struct decoding_context_OA *) decoder->dec_ctx)->sc;
	case BD_DECODER:
		return ((struct decoding_context_BD *) decoder->dec_ctx)->sc;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) decoder->dec_ctx)->sc;
	}
	return 0;
}

int snc_code_overhead(struct snc_decoder *decoder) {
	switch (decoder->d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) decoder->dec_ctx)->overhead;
	case OA_DECODER:
		return ((struct decoding_context_OA *) decoder->dec_ctx)->overhead;
	case BD_DECODER:
		return ((struct decoding_context_BD *) decoder->dec_ctx)->overhead;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) decoder->dec_ctx)->overhead;
	}
	return 0;
}

long long snc_decode_cost(struct snc_decoder *decoder) {
	switch (decoder->d_type) {
	case GG_DECODER:
		return ((struct decoding_context_GG *) decoder->dec_ctx)->operations;
	case OA_DECODER:
		return ((struct decoding_context_OA *) decoder->dec_ctx)->operations;
	case BD_DECODER:
		return ((struct decoding_context_BD *) decoder->dec_ctx)->operations;
	case CBD_DECODER:
		return ((struct decoding_context_CBD *) decoder->dec_ctx)->operations;
	}
	return 0;
}

void snc_free_decoder(struct snc_decoder *decoder) {
	switch (decoder->d_type) {
	case GG_DECODER:
        free_dec_context_GG(((struct decoding_context_GG *) decoder->dec_ctx));
		break;
	case OA_DECODER:
        free_dec_context_OA(((struct decoding_context_OA *) decoder->dec_ctx));
		break;
	case BD_DECODER:
        free_dec_context_BD(((struct decoding_context_BD *) decoder->dec_ctx));
		break;
	case CBD_DECODER:
        free_dec_context_CBD(((struct decoding_context_CBD *) decoder->dec_ctx));
		break;
	}
	decoder->dec_ctx = NULL;
    free(decoder);
    decoder = NULL;
	return;
}
