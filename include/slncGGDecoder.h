#ifndef _CCFD_DECODER_H
#define _CCFD_DECODER_H
#include "slncEncoder.h"

#define FB_THOLD	1
typedef struct node      ID;
typedef struct node_list ID_list;

typedef unsigned long FLAGS;					/* bit flags */
#define _BIT_MASK(n) ( 1UL << (n) )		/* only bit at n is on */
#define _FLAG_SET(x, n) (x) |= _BIT_MASK(n);
#define _FLAG_ON(x, n)  (((x) & _BIT_MASK(n)) == _BIT_MASK(n))
#define _FLAG_OFF(x, n) (((x) & _BIT_MASK(n)) == 0)

struct decoding_context_GG {
    struct slnc_context	*sc;						// The file information
    /********************************
     * Used in decoding LDPC pre-code
     ********************************/
    GF_ELEMENT **evolving_checks;					// Evolving check packets during precode decoding
    // We cannot change decoded packets in file_context directly during iterative decoding,
    // because those packets may be used in decoding other generations later
    int		*check_degrees;							// The number of unknown(undeocded) source neighbors of each check packet

    /******************************
     * Used in decoding generations
     ******************************/
    int finished;
    int decoded;									// record how many packets have been decoded
    int originals;									// record how many source packets are decoded
    struct running_matrix **Matrices;				// record running matrices of each class
    ID_list *recent;								// record most recently decoded packets IDs
    /*******************************************
     * Used if feedback to encoder is allowed
     *******************************************/
    int grecent[FB_THOLD];							// gids of recent FB_THOLD decoded generations (it would wrap around)
    int newgpos;									// pos pointer in grecent where newly decoded gid is stored
    int grcount;									// the number of decoded generations counted from last feedback
    /*******************************************
     * Record for anaylzing decoding performance
     * Not actually needed in implementation
     *******************************************/
    long long operations;							// record the number of computations used
    int overhead;									// record how many packets have been received
};

// to store matrices in processing (needed by the decoder)
struct running_matrix {
    int remaining_rows;								// record how many linearly independent encoding vectors we have
    int remaining_cols;								// record how many source packets remain unknown
    FLAGS erased;									// bits indicating erased columns due to back-substitution 
    GF_ELEMENT **coefficient;
    GF_ELEMENT **message;
};

void create_decoding_context_GG(struct decoding_context_GG *dec_ctx, long datasize, struct slnc_parameter sp);
void free_decoding_context_GG(struct decoding_context_GG *dec_ctx);
void process_packet_GG(struct decoding_context_GG *dec_ctx, struct coded_packet *pkt);
#endif /* _CCFD_DECODER_H */
