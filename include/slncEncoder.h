#ifndef SLNC_ENCODER
#define SLNC_ENCODER
#include <stdio.h>
#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
struct node;
struct node_list;
struct bipartite_graph;

/*
 * Type of SLNC codes
 * 	RAND     - Packets are pesudo-randomly grouped
 * 	BAND     - Packets are grouped to be consecutively overlapped
 * 	WINDWRAP - Similar as BAND, but has wrap around in encoding vectors
 */
#define RAND_SLNC	    0
#define	BAND_SLNC	    1
#define WINDWRAP_SLNC  2

struct source_packet {
    int			id;
    GF_ELEMENT	*syms;			// SIZE_P source message symbols
};

struct slnc_parameter {
    double pcrate;
    int	   size_b;
    int    size_g;
    int    size_p;
    int    type;
};

// Metainfo of the data to be slnc-coded
struct slnc_metainfo {
    long	datasize;					// True data size in bytes. Encoding may append zero packets for alignment.
    double  pcrate;						// precode rate
    int     size_b;
    int		size_g;
    int		size_p;
    int     type;						// Code type
    int		snum;						// Number of source packets splitted from the data.
    int		cnum;						// Number of parity-check packets appended to source packets
    int		gnum;						// Number of subgenerations
};

// Encoding context
struct slnc_context {
    struct  slnc_metainfo   meta;
    struct  subgeneration   **gene; 	// array of pointers each points to a subgeneration.	
    struct  bipartite_graph *graph;
    GF_ELEMENT              **pp;		// Pointers to precoded source packets
};

struct slnc_packet {
    int 		gid;					// subgeneration id;
    GF_ELEMENT	*coes;					// SIZE_G coding coefficients of coded packet
    GF_ELEMENT	*syms;					// SIZE_P symbols of coded packet
};

/**
 * Source packets are grouped into subsets, referred to as 
 * subgenerations in this library to emphasize the 'sub-' concept
 * in sparse code. (Subsets are also referred to as batches, classes, 
 * chunks, segments, or generations in many other network coding literature).
 **/
struct subgeneration {
    int gid;
    int *pktid;							// SIZE_G source packet IDs
};

int slnc_create_enc_context(char *buf, long datasize, struct slnc_context **sc, struct slnc_parameter sp);
int slnc_create_enc_context_from_file(FILE *fp, struct slnc_context **sc, struct slnc_parameter sp);
int slnc_load_file_to_context(FILE *fp, struct slnc_context *sc);
int slnc_free_enc_context(struct slnc_context *sc);
unsigned char *slnc_recover_data(struct slnc_context *sc);
long slnc_recover_data_to_file(FILE *fp, struct slnc_context *sc);
struct slnc_packet *slnc_alloc_empty_packet(int size_g, int size_p);
struct slnc_packet *slnc_generate_packet(struct slnc_context *sc);
int slnc_generate_packet_im(struct slnc_context *sc, struct slnc_packet *pkt);
void slnc_free_packet(struct slnc_packet *pkt);
void print_code_summary(struct slnc_metainfo *meta, int overhead, long long operations);
#endif /* GNC_ENCODER_H */
