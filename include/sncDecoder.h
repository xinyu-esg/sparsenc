#ifndef SNC_DECODER
#define SNC_DECODER
#include "sncEncoder.h"

#define GG_DECODER	0
#define OA_DECODER	1
#define BD_DECODER	2
#define	CBD_DECODER	3

struct snc_decoder;

/**
 * Create an snc decoder given code parameter and decoder type
 *   4 decoders are supported:
 *      GG_DECODER
 *      OA_DECODER
 *      BD_DECODER
 *      CBD_DECODER
 */
struct snc_decoder *snc_create_decoder(struct snc_parameter sp, int d_type);

// Get the encode context that the decoder is working on.
struct snc_context *snc_get_enc_context(struct snc_decoder *decoder, int d_type);

// Feed decoder with snc_packet
void snc_process_packet(struct snc_decoder *decoder, int d_type, struct snc_packet *pkt);

// Check whether the decoder is finished
int snc_decoder_finished(struct snc_decoder *decoder, int d_type); 

// Overhead of code
// Returns the number of received packets of the decoder
int snc_code_overhead(struct snc_decoder *decoder, int d_type);

// Decode cost of the decoder
// Returns the number of finite field operations the decoder has performed
long long snc_decode_cost(struct snc_decoder *decoder, int d_type);

// Free decoder memory
void snc_free_decoder(struct snc_decoder *decoder, int d_type);
#endif
