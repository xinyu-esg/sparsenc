#include "sparsenc.h"
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
    int *ctoo_r;                // record the mapping from current row index to the original row id
    int *ctoo_c;                // record the mapping from current col index to the original row id

    /*performance index*/
    int overhead;               // record how many packets have been received
    int *overheads;             // record how many packets have been received
    long long operations;       // record the number of computations used
};

struct decoding_context_BD *create_dec_context_BD(struct snc_parameters *sp);
void process_packet_BD(struct decoding_context_BD *dec_ctx, struct snc_packet *pkt);
void free_dec_context_BD(struct decoding_context_BD *dec_ctx);

/**
 * File format to store ongoing decoding context
 *
 * snc_parameter
 * decoder_type
 * decoding_context_BD (excluding snc_context)
 *
 */
long save_dec_context_BD(struct decoding_context_BD *dec_ctx, const char *filepath);
struct decoding_context_BD *restore_dec_context_BD(const char *filepath);
