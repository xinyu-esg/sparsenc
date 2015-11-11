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

struct slnc_parameter {
    double pcrate;
    int	   size_b;
    int    size_g;
    int    size_p;
    int    type;
};

// Metainfo of the data to be slnc-coded
struct slnc_metainfo {
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
struct slnc_context {
    struct  slnc_metainfo   meta;
    struct  subgeneration   **gene; 	// array of pointers each points to a generation.	
    struct  bipartite_graph *graph;
    GF_ELEMENT              **pp;		// Pointers to precoded source packets
    //	sp[i][j] - j-th symbol of the i-th source packet
};

struct coded_packet {
    int 		gid;					// generation id;
    GF_ELEMENT	*coes;					// SIZE_G coding coefficients of coded packet
    GF_ELEMENT	*syms;					// SIZE_P symbols of coded packet
};

struct subgeneration {
    int gid;
    int *pktid;							// SIZE_G source packet IDs
};

int create_slnc_context(char *buf, long datasize, struct slnc_context **sc, struct slnc_parameter sp);
int create_slnc_context_from_file(FILE *fp, struct slnc_context **sc, struct slnc_parameter sp);
int load_file_to_slnc_context(FILE *fp, struct slnc_context *sc);
int free_slnc_context(struct slnc_context *sc);
unsigned char *recover_data(struct slnc_context *sc);
long recover_data_to_file(FILE *fp, struct slnc_context *sc);
struct coded_packet *alloc_empty_packet(int size_g, int size_p);
struct coded_packet *generate_slnc_packet(struct slnc_context *sc);
int generate_slnc_packet_im(struct slnc_context *sc, struct coded_packet *pkt);
void free_slnc_packet(struct coded_packet *pkt);
void print_code_summary(struct slnc_metainfo *meta, int overhead, long long operations);
#endif /* GNC_ENCODER_H */
