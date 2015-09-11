/*
 * This file contains functions for erasure-correction coding.
 */
#include "common.h"
#include "gncEncoder.h"
#include <math.h>
#include <sys/stat.h>
static int create_context_from_meta(struct gnc_context *gc);
static int verify_code_parameter(struct gnc_metainfo *meta);
static void perform_precoding(struct gnc_context *gc);
static int group_packets_rand(struct gnc_context *gc);
static int group_packets_band(struct gnc_context *gc);
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
int create_gnc_context(char *buf, long datasize, struct gnc_context **gc, int s_b, int s_g, int s_p, int type)
{
	static char fname[] = "create_gnc_context";
	// Allocate file_context
	if ( (*gc = calloc(1, sizeof(struct gnc_context))) == NULL ) {
		fprintf(stderr, "%s: calloc file_context\n", fname);
		return(-1);
	}
	
	(*gc)->meta.size_b = s_b;
	(*gc)->meta.size_g = s_g;
	(*gc)->meta.size_p = s_p;
	(*gc)->meta.type   = type;

	// Determine packet and generation numbers
	int num_src = ALIGN(datasize, (*gc)->meta.size_p);
	int num_chk = number_of_checks(num_src);
	(*gc)->meta.datasize = datasize;
	(*gc)->meta.snum  = num_src;							  // Number of source packets
	(*gc)->meta.cnum  = num_chk;							  // Number of check packets
	(*gc)->meta.gnum  = ALIGN( (num_src+num_chk), (*gc)->meta.size_b); // Number of disjoint generations grouped from packets
	
	/*
	 * Verify code parameter
	 */
	if (verify_code_parameter(&((*gc)->meta)) != 0) {
		fprintf(stderr, "%s: code parameter is invalid.\n", fname);
		return(-1);
	}
	/*
	 * Create generations, bipartite graph
	 */
	if (create_context_from_meta(*gc) != 0) {
		fprintf(stderr, "%s: create_context_from_meta\n", fname);
		return(-1);
	}

	// Allocating pointers to data
	if (((*gc)->pp = calloc((*gc)->meta.snum+(*gc)->meta.cnum, sizeof(GF_ELEMENT*))) == NULL) {
		fprintf(stderr, "%s: calloc (*gc)->pp\n", fname);
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


int create_gnc_context_from_file(FILE *fp, struct gnc_context **gc, int s_b, int s_g, int s_p, int type)
{
	static char fname[] = "create_gnc_context_from_file";
	// Allocate file_context
	if ( (*gc = calloc(1, sizeof(struct gnc_context))) == NULL ) {
		fprintf(stderr, "%s: calloc file_context\n", fname);
		return(-1);
	}
	
	(*gc)->meta.size_b = s_b;
	(*gc)->meta.size_g = s_g;
	(*gc)->meta.size_p = s_p;
	(*gc)->meta.type   = type;

	// Get file size
	/* Seek to file end */
	if (fseek(fp, 0, SEEK_END) == -1) 
		fprintf(stderr, "%s: fseek SEEK_END\n", fname);
	long datasize = ftell(fp);			/* get file size */
	/* Seek back to file start */
	if (fseek(fp, 0, SEEK_SET) == -1) 
		fprintf(stderr, "%s: fseek SEEK_SET\n", fname);
    /* determine packet and generation numbers */
	int num_src = ALIGN(datasize, (*gc)->meta.size_p);
	int num_chk = number_of_checks(num_src);
	(*gc)->meta.datasize = datasize;
	(*gc)->meta.snum  = num_src;							  // Number of source packets
	(*gc)->meta.cnum  = num_chk;							  // Number of check packets
	(*gc)->meta.gnum  = ALIGN( (num_src+num_chk), (*gc)->meta.size_b); // Number of disjoint generations grouped from packets
	/*
	 * Verify code parameter
	 */
	if (verify_code_parameter(&((*gc)->meta)) != 0) {
		fprintf(stderr, "%s: code parameter is invalid.\n", fname);
		return(-1);
	}
	/*
	 * Create generations, bipartite graph
	 */
	if (create_context_from_meta(*gc) != 0) {
		fprintf(stderr, "%s: create_context_from_meta failed\n", fname);
		return(-1);
	}

	// Allocating pointers to data
	if (((*gc)->pp = calloc((*gc)->meta.snum+(*gc)->meta.cnum, sizeof(GF_ELEMENT*))) == NULL) {
		fprintf(stderr, "%s: calloc (*gc)->pp\n", fname);
		return(-1);
	}

	// Creating context with data read from FILE *fp
	int alread = 0;
	for (int i=0; i<(*gc)->meta.snum+(*gc)->meta.cnum; i++) {
		(*gc)->pp[i] = calloc((*gc)->meta.size_p, sizeof(GF_ELEMENT));
        int toread = (alread+(*gc)->meta.size_p) <= (*gc)->meta.datasize ? (*gc)->meta.size_p : (*gc)->meta.datasize-alread;
		if (fread((*gc)->pp[i], sizeof(GF_ELEMENT), toread, fp) != toread)
			fprintf(stderr, "%s: fread (*gc)->pp[%d]\n", fname, i);
        alread += toread;
    }
	perform_precoding(*gc);

	return(0);
}

static int verify_code_parameter(struct gnc_metainfo *meta)
{
	if (meta->size_b > meta->size_g) {
		fprintf(stderr, "code parameter error: size_b > size_g\n");
		return(-1);
	}
	if (meta->size_b*meta->size_p > meta->datasize) {
		fprintf(stderr, "code parameter error: size_b X size_p > datasize\n");
		return(-1);
	}
	return(0);
}

/*
 * Create gnc context using metadata in fc
 */
static int create_context_from_meta(struct gnc_context *gc)
{
	static char fname[] = "create_gnc_context_meta";
	// Inintialize generation structures
	gc->gene  = malloc(sizeof(struct generation *) * gc->meta.gnum);
	if ( gc->gene == NULL ) {
		fprintf(stderr, "%s: malloc gc->gene\n", fname);
		return(-1);
	}
	for (int j=0; j<gc->meta.gnum; j++) {
		gc->gene[j] = malloc(sizeof(struct generation));
		if ( gc->gene[j] == NULL ) { 
			fprintf(stderr, "%s: malloc gc->gene[%d]\n", fname, j);
			return(-1);
		}
		gc->gene[j]->gid = -1;
		gc->gene[j]->pktid = malloc(sizeof(int)*gc->meta.size_g);
		if ( gc->gene[j]->pktid == NULL ) {
			fprintf(stderr, "%s: malloc gc->gene[%d]->pktid\n", fname, j);
			return(-1);
		}
		memset(gc->gene[j]->pktid, -1, sizeof(int)*gc->meta.size_g);
	}

	int coverage;
	if (gc->meta.type == RAND_GNC_CODE) 
		coverage = group_packets_rand(gc);
	else if (gc->meta.type == BAND_GNC_CODE)
		coverage = group_packets_band(gc);

#if defined(GNCTRACE)
	printf("Data Size: %ld\t Source Packets: %d\t Check Packets: %d\t Generations: %d\t Coverage: %d\n",gc->meta.datasize, gc->meta.snum, gc->meta.cnum,	gc->meta.gnum, coverage);
#endif
	// Creating bipartite graph of the precode
	if (gc->meta.cnum != 0) {
		if ( (gc->graph = malloc(sizeof(BP_graph))) == NULL ) {
			fprintf(stderr, "%s: malloc BP_graph\n", fname);
			return (-1);
		}
		create_bipartite_graph(gc->graph, gc->meta.snum, gc->meta.cnum);
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
		fprintf(stderr, "%s: malloc(datasize) failed.\n", fname);
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

/* recover data to file */
long recover_data_to_file(FILE *fp, struct gnc_context *gc)
{
	static char fname[] = "recover_data";
	long datasize = gc->meta.datasize;
	long alwrote = 0;
	long towrite = datasize;
	
#if defined(GNCTRACE)
	printf("Writing to decoded file.\n");
#endif

	int pc = 0;
	while (alwrote < datasize) {
		towrite = ((alwrote + gc->meta.size_p) <= datasize) ? gc->meta.size_p : datasize - alwrote;
		if (fwrite(gc->pp[pc], sizeof(GF_ELEMENT), towrite, fp) != towrite) 
			fprintf(stderr, "%s: fwrite gc->pp[%d]\n", fname, pc);
		pc++;
		alwrote += towrite;
	}
	return alwrote;
}

// perform systematic LDPC precoding against SRC pkt list and results in a LDPC pkt list
static void perform_precoding(struct gnc_context *gc)
{
	static char fname[] = "perform_precoding";
	// Construct Galois field for encoding and decoding
	constructField(GF_POWER);

	int i, j;
	for (i=0; i<gc->meta.cnum; i++) {
		// Encoding check packet according to the LDPC graph
		NBR_node *nb = gc->graph->l_nbrs_of_r[i]->first;
		while(nb != NULL) {
			int sid = nb->data;				// index of source packet
			// XOR information content
			galois_multiply_add_region(gc->pp[i+gc->meta.snum], gc->pp[sid], 1, gc->meta.size_p, GF_POWER);
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
static int group_packets_rand(struct gnc_context *gc)
{
	int num_p = gc->meta.snum + gc->meta.cnum;
	int num_g = gc->meta.gnum;
	
	int *selected = calloc(num_p, sizeof(int));

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

		// fill in the rest of the generation with packets from other generations
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

/*
 * Group packets to generations that overlap head-to-toe. Each generation's
 * encoding coefficients form a band in GDM.
 */
static int group_packets_band(struct gnc_context *gc)
{
	int num_p = gc->meta.snum + gc->meta.cnum;
	int num_g = gc->meta.gnum;
	
	int *selected = calloc(num_p, sizeof(int));

	int i, j;
	int index;
	int leading_pivot = 0;
	for (i=0; i<num_g; i++) {
		gc->gene[i]->gid = i;
		leading_pivot = i * gc->meta.size_b;
		if (leading_pivot > num_p - gc->meta.size_g)
			leading_pivot = num_p - gc->meta.size_g;
		for (j=0; j<gc->meta.size_g; j++) {
			index = leading_pivot + j;
			selected[index] += 1;
			gc->gene[i]->pktid[j] = index;
		}	
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
		co = (GF_ELEMENT) rand() % (1 << GF_POWER);					// Randomly generated coding coefficient
		galois_multiply_add_region(pkt->syms, gc->pp[pktid], co, gc->meta.size_p, GF_POWER);
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
