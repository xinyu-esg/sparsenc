#ifndef PP_DECODER_H
#define PP_DECODER_H
#include "sparsenc.h"

#define FORWARD         1
#define FINALFORWARD    2
#define FINALBACKWARD   3
/*
 * Perpetual code DECODING CONTEXT
 */
struct decoding_context_PP
{
    // GNC context
    struct snc_context *sc;

    int finished;               // an indicator tracking the finish of decoding
    int stage;                  // stage of the decoding
    int pivots;                 // number of pivots that have been received

    // decoding matrix
    struct row_vector **row;    // NUM_PP rows for storing coefficient vectors
    // row[i] represents the i-th row starting from the diagonal element A[i][i]
    GF_ELEMENT **message;       // NUM_PP rows for storing message symbols

    /*performance index*/
    int overhead;               // record how many packets have been received
    long long operations;       // record the number of computations used
};

struct decoding_context_PP *create_dec_context_PP(struct snc_parameters *sp);
void process_packet_PP(struct decoding_context_PP *dec_ctx, struct snc_packet *pkt);
void free_dec_context_PP(struct decoding_context_PP *dec_ctx);

/**
 * File format to store ongoing decoding context
 * abc.decoderCBD.part
 *
 * snc_parameter
 * decoder_type
 * decoding_context_CBD (excluding snc_context)
 *
 */
long save_dec_context_PP(struct decoding_context_PP *dec_ctx, const char *filepath);
struct decoding_context_PP *restore_dec_context_PP(const char *filepath);
#endif
