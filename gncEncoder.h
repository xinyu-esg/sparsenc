#ifndef _GNC_ENCODER_H
#define _GNC_ENCODER_H
#include "common.h"
#include "galois.h"
#include "bipartite.h"

struct source_packet {
	int			id;
	GF_ELEMENT	*syms;					// SIZE_P source message symbols
};

// Metainfo of the data to be gnc-coded
struct gnc_metainfo {
	long	datasize;					// True data size in bytes. GNC may need to append  
										// zero packets to align generations.
	int     size_b;
	int		size_g;
	int		size_p;
	int		snum;						// Number of source packets splitted from the data.
	int		cnum;						// Number of parity-check packets appended to source packets
	int		gnum;						// Number of generations
};

// Encoding context
struct gnc_context {
	struct  gnc_metainfo    meta;
	struct  generation      **gene; 	// array of pointers each points to a generation.	
	struct  bipartite_graph *graph;
	GF_ELEMENT              **pp;		// Pointers to precoded source packets
										//	sp[i][j] - j-th symbol of the i-th source packet
};

struct coded_packet {
	int 		gid;					// generation id;
	GF_ELEMENT	*coes;					// SIZE_G coding coefficients of coded packet
	GF_ELEMENT	*syms;					// SIZE_P symbols of coded packet
};

struct generation {
	int gid;
	int *pktid;							// SIZE_G source packet IDs
};

int create_gnc_context(char *buf, long datasize, struct gnc_context **gc, int s_b, int s_g, int s_p);
int free_gnc_context(struct gnc_context *gc);
unsigned char *recover_data(struct gnc_context *gc);
struct coded_packet *generate_gnc_packet(struct gnc_context *gc);
void free_gnc_packet(struct coded_packet *pkt);
#endif /* _GNC_ENCODER_H */
