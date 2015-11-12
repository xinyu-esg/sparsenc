#include "slncEncoder.h"
/*
 * BD (band GNC code) DECODING CONTEXT
 */
struct slnc_dec_context_BD
{
    // GNC context
    struct slnc_context *sc;

    int finished;				// an indicator tracking the finish of decoding
    int DoF;					// total true DoF that the receiver has received
    int de_precode;				// apply precode or not
    int inactivated;			// total number of inactivated packets among overlapping packets

    // decoding matrix
    GF_ELEMENT **coefficient;						//[NUM_PP][NUM_PP];
    GF_ELEMENT **message;							//[NUM_PP][EXT_N];

    // the following two mappings are to record pivoting processings
    int *otoc_mapping;								//[NUM_PP] record the mapping from original packet id to current column index
    int *ctoo_mapping;								//[NUM_PP] record the mapping from current column index to the original packet id

    /*performance index*/
    int overhead;									// record how many packets have been received
    int *overheads;									// record how many packets have been received
    long long operations;							// record the number of computations used
};

void slnc_create_dec_context_BD(struct slnc_dec_context_BD *dec_ctx, long datasize, struct slnc_parameter sp);
void slnc_process_packet_BD(struct slnc_dec_context_BD *dec_ctx, struct slnc_packet *pkt);
void slnc_free_dec_context_BD(struct slnc_dec_context_BD *dec_ctx);

