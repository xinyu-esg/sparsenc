#ifndef GNC_ENCODER_H
#define GNC_ENCODER_H
#include <stdio.h>
#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
struct node;
struct node_list;
struct bipartite_graph;

/*
 * Type of GNC code
 * 	RAND - Generations are pesudo-randomly grouped
 * 	BAND - Generations are overlapped consecutively
 * 	WNND - Windowed, generations are overlapped consecutively and wrap around
 */
#define RAND_GNC_CODE	    0
#define	BAND_GNC_CODE	    1
#define WINDWRAP_GNC_CODE   2

struct source_packet {
    int			id;
    GF_ELEMENT	*syms;					// SIZE_P source message symbols
};

struct gnc_parameter {
    double pcrate;
    int	   size_b;
    int    size_g;
    int    size_p;
    int    type;
};

// Metainfo of the data to be gnc-coded
struct gnc_metainfo {
    long	datasize;					// True data size in bytes. GNC may append zero packets for alignment.
    double  pcrate;						// precode rate
    int     size_b;
    int		size_g;
    int		size_p;
    int     type;						// GNC code type - how generations are grouped
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

int create_gnc_context(char *buf, long datasize, struct gnc_context **gc, struct gnc_parameter gp);
int create_gnc_context_from_file(FILE *fp, struct gnc_context **gc, struct gnc_parameter gp);
int load_file_to_gnc_context(FILE *fp, struct gnc_context *gc);
int free_gnc_context(struct gnc_context *gc);
unsigned char *recover_data(struct gnc_context *gc);
long recover_data_to_file(FILE *fp, struct gnc_context *gc);
struct coded_packet *alloc_empty_packet(int size_g, int size_p);
struct coded_packet *generate_gnc_packet(struct gnc_context *gc);
int generate_gnc_packet_im(struct gnc_context *gc, struct coded_packet *pkt);
void free_gnc_packet(struct coded_packet *pkt);
void print_code_summary(struct gnc_metainfo *meta, int overhead, long long operations);
#endif /* GNC_ENCODER_H */
