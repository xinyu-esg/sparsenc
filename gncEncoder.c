/*
 * This file contains functions for erasure-correction coding.
 */
#include "common.h"
#include "gncEncoder.h"
#include <math.h>
#include <sys/stat.h>
static void perform_precoding(struct gnc_context *gc);
static int group_packets(struct gnc_context *gc);
static void encode_packet(struct gnc_context *gc, int gid, struct coded_packet *pkt);
static int schedule_generation(struct gnc_context *gc);
static int is_prime(int number); 
static int number_of_checks(int snum); 
/*
 * Create a GNC context containing meta information about the data to be encoded.
 *   buf      - Buffer containing bytes of data to be encoded
 *   datasize - size of data in bytes to be encoded
 *   gc       - pointer to the gnc_context where the context will be stored
 *   s_b	  - base generation size, size_b
 *   s_g	  - full generation size, size_g
 *   s_p	  - number of information symbols, size_p
 *
 * Return
 *   0  - Create successfully
 *   -1 - Create failed
 */
int create_gnc_context(char *buf, long datasize, struct gnc_context **gc, int s_b, int s_g, int s_p)
{
	static char fname[] = "create_gnc_context";
	// Allocate file_context
	if ( (*gc = calloc(1, sizeof(struct gnc_context))) == NULL ) {
		printf("%s: calloc file_context\n", fname);
		return(-1);
	}
	
	(*gc)->meta.size_b = s_b;
	(*gc)->meta.size_g = s_g;
	(*gc)->meta.size_p = s_p;

	// Determine packet and generation numbers
	int num_src = ALIGN(datasize, (*gc)->meta.size_p);
	int num_chk = number_of_checks(num_src);
	(*gc)->meta.datasize = datasize;
	(*gc)->meta.snum  = num_src;							  // Number of source packets
	(*gc)->meta.cnum  = num_chk;							  // Number of check packets
	(*gc)->meta.gnum  = ALIGN( (num_src+num_chk), (*gc)->meta.size_b); // Number of disjoint generations grouped from packets
	
	// Inintialize generation structures
	(*gc)->gene  = malloc(sizeof(struct generation *) * (*gc)->meta.gnum);
	if ( (*gc)->gene == NULL ) {
		printf("%s: malloc (*gc)->gene\n", fname);
		return(-1);
	}
	for (int j=0; j<(*gc)->meta.gnum; j++) {
		(*gc)->gene[j] = malloc(sizeof(struct generation));
		if ( (*gc)->gene[j] == NULL ) { 
			printf("%s: malloc (*gc)->gene[%d]\n", fname, j);
			return(-1);
		}
		(*gc)->gene[j]->gid = -1;
		(*gc)->gene[j]->pktid = malloc(sizeof(int)*(*gc)->meta.size_g);
		memset((*gc)->gene[j]->pktid, -1, sizeof(int)*(*gc)->meta.size_g);
	}
	int coverage = group_packets(*gc);
	printf("Data Size: %ld\t Source Packets: %d\t Check Packets: %d\t Generations: %d\t Coverage: %d\n",(*gc)->meta.datasize, (*gc)->meta.snum, (*gc)->meta.cnum,	(*gc)->meta.gnum, coverage);
	
	// Creating bipartite graph of the precode
	if ((*gc)->meta.cnum != 0) {
		if ( ((*gc)->graph = malloc(sizeof(BP_graph))) == NULL ) {
			printf("%s: malloc BP_graph\n", fname);
			return (-1);
		}
		create_bipartite_graph((*gc)->graph, (*gc)->meta.snum, (*gc)->meta.cnum);
	}

	// Allocating pointers to data
	if (((*gc)->pp = calloc((*gc)->meta.snum+(*gc)->meta.cnum, sizeof(GF_ELEMENT*))) == NULL) {
		printf("%s: calloc (*gc)->pp\n", fname);
		return(-1);
	}

	// Creating context with to-be-encoded data
	if (buf != NULL) {
		int alread = 0;
		for (int i=0; i<(*gc)->meta.snum+(*gc)->meta.cnum; i++) {
        	(*gc)->pp[i] = calloc((*gc)->meta.size_p, sizeof(GF_ELEMENT));
        	int toread = (alread+(*gc)->meta.size_p) <= (*gc)->meta.datasize ? (*gc)->meta.size_p : (*gc)->meta.datasize-alread;
        	memcpy((*gc)->pp[i], buf+alread, toread*sizeof(GF_ELEMENT));
        	alread += toread;
    	}
		perform_precoding(*gc);
	}

	return(0);
}

int free_gnc_context(struct gnc_context *gc)
{
	int i;
	for (i=gc->meta.snum+gc->meta.cnum-1; i>=0; i--) {
		if (gc->pp[i] != NULL) {
			free(gc->pp[i]);
			gc->pp[i] = NULL;
		}
	}
	free(gc->pp);
	for (i=gc->meta.gnum-1; i>=0; i--) {
		free(gc->gene[i]->pktid);			// free packet IDs
		free(gc->gene[i]);					// free generation itself
		gc->gene[i] = NULL;
	}
	free(gc->gene);
	free_bipartite_graph(gc->graph);
	free(gc);
	return(0);
}

unsigned char *recover_data(struct gnc_context *gc)
{
	static char fname[] = "recover_data";
	long datasize = gc->meta.datasize;
	long alwrote = 0;
	long towrite = datasize;
	
	unsigned char *data;
	if ( (data = malloc(datasize)) == NULL) {
		printf("%s: malloc(datasize) failed.\n", fname);
		return NULL;
	}
	int pc = 0;
	while (alwrote < datasize) {
		towrite = ((alwrote + gc->meta.size_p) <= datasize) ? gc->meta.size_p : datasize - alwrote;
		memcpy(data+alwrote, gc->pp[pc++], sizeof(GF_ELEMENT)*towrite);
		alwrote += towrite;
	}
	return data;
}

// perform systematic LDPC precoding against SRC pkt list and results in a LDPC pkt list
static void perform_precoding(struct gnc_context *gc)
{
	static char fname[] = "perform_precoding";
	// Construct Galois field for encoding and decoding
	constructField(GF_ORDER);

	int i, j;
	for (i=0; i<gc->meta.cnum; i++) {
		// Encoding check packet according to the LDPC graph
		NBR_node *nb = gc->graph->l_nbrs_of_r[i]->first;
		while(nb != NULL) {
			int sid = nb->data;				// index of source packet
			// XOR information content
			galois_multiply_add_region(gc->pp[i+gc->meta.snum], gc->pp[sid], 1, gc->meta.size_p, GF_ORDER);
			// move to next possible neighbour node of current check
			nb = nb->next;
		}
	}
}

/*
 * This routine uses a deterministic grouping scheme, so the need of sending
 * grouping information to clients is removed. The only information clients 
 * need to know is the number of packets, base size, and generation size. 
 */
static int group_packets(struct gnc_context *gc)
{
	int num_p = gc->meta.snum + gc->meta.cnum;
	int num_g = gc->meta.gnum;
	
	int *selected = malloc(sizeof(int)*num_p);
	memset(selected, 0, sizeof(int)*num_p);

	int i, j;
	int index;
	int rotate = 0;
	for (i=0; i<num_g; i++) {
		gc->gene[i]->gid = i;
		// split packets into disjoint groups
		for (j=0; j<gc->meta.size_b; j++) {
			index = (i * gc->meta.size_b + j) % num_p;				// source packet index

			while (has_item(gc->gene[i]->pktid, index, j) != -1)
				index++;
			gc->gene[i]->pktid[j] = index;
			selected[index] += 1;
		}

		for (j=gc->meta.size_b; j<gc->meta.size_g; j++) {
			int tmp = i - (j - gc->meta.size_b + 7);
			int start = tmp >= 0 ? tmp : tmp+num_g;
			if (start == i)
				start++;
			index = (start * gc->meta.size_b + (j - gc->meta.size_b + rotate) % (gc->meta.size_g)) % num_p;
			while (has_item(gc->gene[i]->pktid, index, j) != -1)
				index++;
			gc->gene[i]->pktid[j] = index;
			selected[index] += 1;
		}
		rotate = (rotate + 7) % (gc->meta.size_g);
	}
	int coverage = 0;
	for (i=0; i<num_p; i++)
		coverage += selected[i];

	free(selected);
	return coverage;
}

struct coded_packet *generate_gnc_packet(struct gnc_context *gc)
{
	struct coded_packet *pkt = malloc(sizeof(struct coded_packet));
	if (pkt == NULL)
		return NULL;
	pkt->coes = calloc(gc->meta.size_g, sizeof(GF_ELEMENT));
	if (pkt->coes == NULL) {
		free(pkt);
		return NULL;
	}
	pkt->syms = calloc(gc->meta.size_p, sizeof(GF_ELEMENT));
	if (pkt->syms == NULL) {
		free(pkt);
		return NULL;
	}
	int gid = schedule_generation(gc);
	encode_packet(gc, gid, pkt);
	return pkt;
}

void free_gnc_packet(struct coded_packet *pkt)
{
	free(pkt->coes);
	free(pkt->syms);
	free(pkt);
}


static void encode_packet(struct gnc_context *gc, int gid, struct coded_packet *pkt)
{
	pkt->gid = gid;
	int i;
	GF_ELEMENT co;
	int pktid;
	for (i=0; i<gc->meta.size_g; i++) {
		pktid = gc->gene[gid]->pktid[i];							// The i-th packet of the gid-th generation
		co = (GF_ELEMENT) rand() % (1 << GF_ORDER);					// Randomly generated coding coefficient
		galois_multiply_add_region(pkt->syms, gc->pp[pktid], co, gc->meta.size_p, GF_ORDER);
		pkt->coes[i] = co;
	}
}

static int schedule_generation(struct gnc_context *gc)
{
	int gid = rand() % (gc->meta.gnum);
	return gid;
}

static int is_prime(int number) 
{
    int i;
    for (i=2; i*i<=number; i++) {
        if (number % i == 0)
			return 0;
    }
    return 1;
}

// Return the number of required LDPC check symbols given the number of source packets.
static int number_of_checks(int snum) 
{
	int x = (int) floor( sqrt( 2 * snum ) );
	while ( x * (x - 1) < 2 * snum )
		x++;

	int c = (int) ceil( 0.01 * snum ) + x;
	while ( !is_prime(c++) )
		;

	return c;
}