#ifndef CBD_DECODER_H
#define CBD_DECODER_H
#include "sparsenc.h"

/* Row vector of a matrix */
struct row_vector
{
    int len;            // length of the row
    GF_ELEMENT *elem;   // elements of the row
};
/*
 * Compact BD (band GNC code) DECODING CONTEXT
 */
struct decoding_context_CBD
{
    // GNC context
    struct snc_context *sc;

    int finished;               // an indicator tracking the finish of decoding
    int DoF;                    // total true DoF that the receiver has received
    int de_precode;             // apply precode or not
    int naive;                  // decode in naive mode (for non-band code)

    // decoding matrix
    struct row_vector **row;    // NUM_PP rows for storing coefficient vectors
    // row[i] represents the i-th row starting from the diagonal element A[i][i]
    GF_ELEMENT **message;       // NUM_PP rows for storing message symbols

    /*performance index*/
    int overhead;               // record how many packets have been received
    long long operations;       // record the number of computations used
#if defined(SIMULATION)
    long long ops1;             // operations of forward sub
    long long ops2;             // operations of applying precode
    long long ops3;             // operations of backward sub
#endif
};

struct decoding_context_CBD *create_dec_context_CBD(struct snc_parameters *sp);
void process_packet_CBD(struct decoding_context_CBD *dec_ctx, struct snc_packet *pkt);
void free_dec_context_CBD(struct decoding_context_CBD *dec_ctx);

/**
 * File format to store ongoing decoding context
 * abc.decoderCBD.part
 *
 * snc_parameter
 * decoder_type
 * decoding_context_CBD (excluding snc_context)
 *
 */
long save_dec_context_CBD(struct decoding_context_CBD *dec_ctx, const char *filepath);
struct decoding_context_CBD *restore_dec_context_CBD(const char *filepath);
#endif
