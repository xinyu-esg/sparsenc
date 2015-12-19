#include "snc.h"
/*
 * BD (band GNC code) DECODING CONTEXT
 */
struct decoding_context_BD
{
    // GNC context
    struct snc_context *sc;

    int finished;               // an indicator tracking the finish of decoding
    int DoF;                    // total true DoF that the receiver has received
    int de_precode;             // apply precode or not
    int inactivated;            // total number of inactivated packets among overlapping packets

    // decoding matrix
    GF_ELEMENT **coefficient;   //[NUM_PP][NUM_PP];
    GF_ELEMENT **message;       //[NUM_PP][EXT_N];

    // the following two mappings are to record pivoting processings
    int *otoc_mapping;          //[NUM_PP] record the mapping from original packet id to current column index
    int *ctoo_mapping;          //[NUM_PP] record the mapping from current column index to the original packet id

    /*performance index*/
    int overhead;               // record how many packets have been received
    int *overheads;             // record how many packets have been received
    long long operations;       // record the number of computations used
};

void create_dec_context_BD(struct decoding_context_BD *dec_ctx, struct snc_parameter sp);
void process_packet_BD(struct decoding_context_BD *dec_ctx, struct snc_packet *pkt);
void free_dec_context_BD(struct decoding_context_BD *dec_ctx);

