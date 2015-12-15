#ifndef CBD_DECODER_H
#define CBD_DECODER_H
#include "snc.h"

/* Row vector of a matrix */
struct row_vector
{
    int len;			// length of the row
    GF_ELEMENT *elem;	// elements of the row
};
/*
 * Compact BD (band GNC code) DECODING CONTEXT
 */
struct decoding_context_CBD
{
    // GNC context
    struct snc_context *sc;

    int finished;				// an indicator tracking the finish of decoding
    int DoF;					// total true DoF that the receiver has received
    int de_precode;				// apply precode or not

    // decoding matrix
    struct row_vector **row;	// NUM_PP rows for storing coefficient vectors
    // row[i] represents the i-th row starting from the diagonal element A[i][i]
    GF_ELEMENT **message;		// NUM_PP rows for storing message symbols

    /*performance index*/
    int overhead;				// record how many packets have been received
    long long operations;		// record the number of computations used
};

void create_dec_context_CBD(struct decoding_context_CBD *dec_ctx, struct snc_parameter sp);
void process_packet_CBD(struct decoding_context_CBD *dec_ctx, struct snc_packet *pkt);
void free_dec_context_CBD(struct decoding_context_CBD *dec_ctx);
#endif
