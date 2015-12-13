#ifndef SNC_ENCODER
#define SNC_ENCODER
#include <stdio.h>
#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
/*
 * Type of SNC codes
 * 	RAND     - Packets are pesudo-randomly grouped
 * 	BAND     - Packets are grouped to be consecutively overlapped
 * 	WINDWRAP - Similar as BAND, but has wrap around in encoding vectors
 */
#define RAND_SNC        0
#define	BAND_SNC        1
#define WINDWRAP_SNC    2

// Sparse network code encode context
struct snc_context;

struct snc_packet {
    int 		gid;					// subgeneration id;
    GF_ELEMENT	*coes;					// SIZE_G coding coefficients of coded packet
    GF_ELEMENT	*syms;					// SIZE_P symbols of coded packet
};

struct snc_parameter {
    long   datasize;
    double pcrate;
    int	   size_b;
    int    size_g;
    int    size_p;
    int    type;
};

// Metainfo of the data to be snc-coded
struct snc_metainfo {
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

/**
 * Create encode context from a message buffer pointed by buf. Code parameters
 * are provided in sp.
 *
 * Return Values:
 *   On success, a pointer to the allocated encode context is returned;
 *   On error, NULL is returned, and errno is set appropriately.
 **/
struct snc_context *snc_create_enc_context(char *buf, struct snc_parameter sp);

struct snc_metainfo *snc_get_metainfo(struct snc_context *sc);

int snc_load_file_to_context(const char *filepath, long start, struct snc_context *sc);

int snc_free_enc_context(struct snc_context *sc);

unsigned char *snc_recover_data(struct snc_context *sc);

long snc_recover_to_file(const char *filepath, struct snc_context *sc);

struct snc_packet *snc_alloc_empty_packet(int size_g, int size_p);

struct snc_packet *snc_generate_packet(struct snc_context *sc);

int snc_generate_packet_im(struct snc_context *sc, struct snc_packet *pkt);

void snc_free_packet(struct snc_packet *pkt);

void print_code_summary(struct snc_context *sc, int overhead, long long operations);
#endif /* SNC_ENCODER */
