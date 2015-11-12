#ifndef CBD_DECODER
#define CBD_DECODER
#include "slncEncoder.h"

/* Row vector of a matrix */
struct row_vector
{
    int len;			// length of the row
    GF_ELEMENT *elem;	// elements of the row
};
/*
 * Compact BD (band GNC code) DECODING CONTEXT
 */
struct slnc_dec_context_CBD
{
    // GNC context
    struct slnc_context *sc;

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

void slnc_create_dec_context_CBD(struct slnc_dec_context_CBD *dec_ctx, long datasize, struct slnc_parameter sp);
void slnc_process_packet_CBD(struct slnc_dec_context_CBD *dec_ctx, struct slnc_packet *pkt);
void slnc_free_dec_context_CBD(struct slnc_dec_context_CBD *dec_ctx);
#endif
