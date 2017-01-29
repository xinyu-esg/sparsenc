/*--------------------- common.h ------------------------
 *  Internal header file
 *------------------------------------------------------*/
#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "sparsenc.h"

/* log levels */
#define TRACE       5

#define ALIGN(a, b) ((a) % (b) == 0 ? (a)/(b) : (a)/(b) + 1)
#define RESIDUAL(a, b) ((b) * ALIGN((a), (b)) - (a))

#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
typedef struct node      NBR_node;
typedef struct node_list NBR_nodes;
// node of singly linked list
struct node {
    int data;
    GF_ELEMENT ce;     // coefficient associated with node's data (for bipartite graph)
    struct node *next;
};

struct node_list {
    struct node *first;
    struct node *last;
};

// Bipartitle graph for LDPC code
typedef struct bipartite_graph {
    int         nleft;
    int         nright;
    int         binaryce;       // Whether coefficients of edges are 1 or higher order
    NBR_nodes **l_nbrs_of_r;    // left side neighbours of right
    NBR_nodes **r_nbrs_of_l;    // right side neighbours of left
} BP_graph;

/**
 * Source packets are grouped into subsets, referred to as
 * subgenerations in this library to emphasize the 'sub-' concept
 * in sparse code. (Subsets are also referred to as batches, classes,
 * chunks, segments, or generations in many other network coding literature).
 **/
struct subgeneration {
    int gid;
    int *pktid;                 // SIZE_G source packet IDs
};

/**
 * Definition of snc_context
 **/
struct snc_context {
    struct  snc_parameters    params;
    int                       snum;     // Number of source packets splitted
    int                       cnum;     // Number of parity-checks(cnum ~= snum * pcrate)
    int                       gnum;     // Number of subgenerations
    struct  subgeneration   **gene;     // array of pointers each points to a subgeneration.
    struct  bipartite_graph  *graph;
    GF_ELEMENT              **pp;       // Pointers to precoded source packets
    int                      *nccount;  // Count of coded packets generated from each subgeneration
    int                       count;    // Count of total coded packets generated
};


/*
 * Buffer for storing SNC packets (for recoding)
 *
 * Buffer size specifies how many packets are saved for
 * each subgeneration. "FIFO" strategy is used when buffer
 * size is reached; the oldest buffered packet will be
 * discarded when a subgeneration buffer is full while a new
 * packet belonging to the subgeneration arrives.
 *
 * Buffer data structure
 *
 * gbuf --> gbuf[0]
 *                     snc_packet  snc_packet ...
 *          gbuf[1]          ^            ^
 *                           |            |
 *          gbuf[2] --> gbuf[2][0]   gbuf[2][1] ....
 *            .
 *            .
 *            .
 *
 * Systematic packet buffer (systematic code)
 *        snc_packet
 *             ^
 *             |
 * sysbuf[0] sysbuf[1] sysbuf[2] ... 
 */
struct snc_buffer {
    struct snc_parameters  params;  // Meta info of the code
    int                    snum;    // Number of source packets
    int                    cnum;    // Number of parity-check packets
    int                    gnum;    // Number of subgenerations
    int                    size;    // Number of bufferred packets of each subgeneration
    int                    nemp;    // Number of non-empty subgeneration buffers
    struct snc_packet   ***gbuf;    // Pointers to subgeneration buffers
    int                   *nc;      // Number of currently buffered packets of each generation
    int                   *pn;      // Positions to store next packet of each subgeneration
    int                   *nsched;  // Number of scheduled times of each subgeneration
    // This is used during systematic scheduling
    //int                   *prevuc;    // position of last scheduled uncoded packet
    //int                   *lastuc;  // position of last buffered uncoded packet
    struct snc_packet    **sysbuf;   // Buffered uncoded packet (needed for systematic code)
    int                    sysnum;  // number of buffered systematic packet
    int                    sysptr;  // pointer of already scheduled systematic packet
};

/* Row vector of a matrix */
struct row_vector
{
    int len;            // length of the row
    GF_ELEMENT *elem;   // elements of the row
};


/* common.c */
void set_loglevel(const char *level);
int get_loglevel();
int has_item(int array[], int item, int length);
void append_to_list(struct node_list *list, struct node *nd);
int remove_from_list(struct node_list *list, int data);
int exist_in_list(struct node_list *list, int data);
void clear_list(struct node_list *list);
void free_list(struct node_list *list);
unsigned char get_bit_in_array(unsigned char *coes, int i);
void set_bit_in_array(unsigned char *coes, int i);
//int snc_rand(void);
//void snc_srand(unsigned int seed);
/* bipartite.c */
int number_of_checks(int snum, double r);
int create_bipartite_graph(BP_graph *graph, int nleft, int nright);
void free_bipartite_graph(BP_graph *graph);
// mt19937ar.c
void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
#endif /* COMMON_H */
