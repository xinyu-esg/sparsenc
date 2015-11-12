#ifndef OA_DECODER
#define OA_DECODER
#include "slncEncoder.h"
/*
 * OA (overlap aware) DECODING CONTEXT is used at each destination node
 * to record the status of the decoder. It is the core
 * information needed by the decoder
 */
struct slnc_dec_context_OA
{
    // GNC context
    struct slnc_context *sc;

    int aoh;										// Allowed number of overhead packets (OHS)
    int finished;									// an indicator tracking the finish of decoding
    int OA_ready;									// the decoder has reached the necessary condition for OA decoding
    int local_DoF;									// total DOF that have been received within generations, total_DoF==NUM_SRC, then decodable
    int global_DoF;									// total true DoF that the receiver has received

    // Local decoding matrices
    struct running_matrix **Matrices;				//[CLASS_NUM] record running matrices of each class
    // Global decoding matrix
    GF_ELEMENT **JMBcoefficient;					//[NUM_SRC+OHS+CHECKS][NUM_PP];
    GF_ELEMENT **JMBmessage;						//[NUM_SRC+OHS+CHECKS][EXT_N];

    // the following two mappings are to record pivoting processings
    int *otoc_mapping;								//[NUM_PP] record the mapping from original packet id to current column index
    int *ctoo_mapping;								//[NUM_PP] record the mapping from current column index to the original packet id
    int inactives;									// total number of inactivated packets among overlapping packets
    //SRC_PKT				*dpkt[NUM_PP];				// decoded packet in the original order

    int overhead;									// record how many packets have been received
    long long operations;							// record the number of computations used
};


// to store matrices in processing (needed by the decoder)
struct running_matrix
{
    GF_ELEMENT **coefficient;						//[CLASS_SIZE][CLASS_SIZE];		
    GF_ELEMENT **message;							//[CLASS_SIZE][EXT_N];
};

void slnc_create_dec_context_OA(struct slnc_dec_context_OA *dec_ctx, long datasize, struct slnc_parameter sp, int aoh);
void slnc_process_packet_OA(struct slnc_dec_context_OA *dec_ctx, struct slnc_packet *pkt);
void slnc_free_dec_context_OA(struct slnc_dec_context_OA *dec_ctx);
#endif
