#ifndef SLNC_DECODER
#define SLNC_DECODER
#include "slncEncoder.h"

#define GG_DECODER	0
#define OA_DECODER	1
#define BD_DECODER	2
#define	CBD_DECODER	3

int slnc_create_dec_context(int d_type, long datasize, struct slnc_parameter sp);
void slnc_process_packet(int d_type, struct slnc_packet *pkt);
int slnc_is_dec_finished(int d_type); 
struct slnc_context *slnc_decoded_context(int d_type);
int slnc_dec_overhead(int d_type);
long long slnc_dec_operations(int d_type);
void slnc_free_dec_context(int d_type);
#endif
