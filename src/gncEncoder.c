/**************************************************************
 * 		gncEncoder.c
 *
 * Functions for GNC encoding. Coded packets can be generated
 * from memory buffer or files.
 **************************************************************/
#include "common.h"
#include "galois.h"
#include "bipartite.h"
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
int create_gnc_context(char *buf, long datasize, struct gnc_context **gc, struct gnc_parameter gp)
{
	static char fname[] = "create_gnc_context";
	// Allocate file_context
	if ( (*gc = calloc(1, sizeof(struct gnc_context))) == NULL ) {
		fprintf(stderr, "%s: calloc file_context\n", fname);
		return(-1);
	}
	(*gc)->meta.pcrate   = gp.pcrate;	
	(*gc)->meta.size_b   = gp.size_b;
	(*gc)->meta.size_g   = gp.size_g;
	(*gc)->meta.size_p   = gp.size_p;
	(*gc)->meta.type     = gp.type;

	// Determine packet and generation numbers
	int num_src = ALIGN(datasize, (*gc)->meta.size_p);
	int num_chk = number_of_checks(num_src, (*gc)->meta.pcrate);
	(*gc)->meta.datasize = datasize;
	(*gc)->meta.snum  = num_src;							  // Number of source packets
	(*gc)->meta.cnum  = num_chk;							  // Number of check packets
	if ((*gc)->meta.type == BAND_GNC_CODE) 
		(*gc)->meta.gnum  = ALIGN((num_src+num_chk-(*gc)->meta.size_g), (*gc)->meta.size_b) + 1; 
	else
		(*gc)->meta.gnum  = ALIGN( (num_src+num_chk), (*gc)->meta.size_b); 
	
	/*
	 * Verify code gpmeter
	 */
	if (verify_code_parameter(&((*gc)->meta)) != 0) {
		fprintf(stderr, "%s: code gpmeter is invalid.\n", fname);
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

	/*--------- Creating context with to-be-encoded data ---------------
	 *
	 * Currently we copy data from user-provided buffer. In the future,
	 * we can make most of pp[i] directly point to user-provided buffer
	 * address. Only the last source packet (due to zero-padding) and 
	 * the cnum parity-check packets need to allocate new memory.
	 *
	 *------------------------------------------------------------------*/
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
	// Construct Galois field for encoding and decoding
	constructField(GF_POWER);

	return(0);
}


int create_gnc_context_from_file(FILE *fp, struct gnc_context **gc, struct gnc_parameter gp)
{
	static char fname[] = "create_gnc_context_from_file";
	// Get file size
	/* Seek to file end */
	if (fseek(fp, 0, SEEK_END) == -1) 
		fprintf(stderr, "%s: fseek SEEK_END\n", fname);
	long datasize = ftell(fp);			/* get file size */
	/* Seek back to file start */
	if (fseek(fp, 0, SEEK_SET) == -1) 
		fprintf(stderr, "%s: fseek SEEK_SET\n", fname);

	/* Create gnc_context without actual data */
	create_gnc_context(NULL, datasize, gc, gp);
	return(0);
}

/*
 * Load file data into gnc_context created for the file. It's the
 * caller's responsibility to ensure that the pair of FILE* and 
 * gnc_context matches.
 */
int load_file_to_gnc_context(FILE *fp, struct gnc_context *gc)
{
	static char fname[] = "load_file_to_gnc_context";
	int alread = 0;
	for (int i=0; i<gc->meta.snum+gc->meta.cnum; i++) {
		gc->pp[i] = calloc(gc->meta.size_p, sizeof(GF_ELEMENT));
        int toread = (alread+gc->meta.size_p) <= gc->meta.datasize ? gc->meta.size_p : gc->meta.datasize-alread;
		if (fread(gc->pp[i], sizeof(GF_ELEMENT), toread, fp) != toread) {
			fprintf(stderr, "%s: fread gc->pp[%d]\n", fname, i);
			return (-1);
		}
        alread += toread;
    }
	perform_precoding(gc);
	return (0);
}

static int verify_code_parameter(struct gnc_metainfo *meta)
{
	if (meta->size_b > meta->size_g) {
		fprintf(stderr, "code gpmeter error: size_b > size_g\n");
		return(-1);
	}
	if (meta->size_b*meta->size_p > meta->datasize) {
		fprintf(stderr, "code gpmeter error: size_b X size_p > datasize\n");
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
	if (gc->graph != NULL)
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
		if (leading_pivot > num_p - gc->meta.size_g) {
#if defined(GNCTRACE)
			printf("Band lead of gid: %d is modified\n", i);
#endif
			leading_pivot = num_p - gc->meta.size_g;
		}
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

/*
 * Allocate an empty GNC coded packet
 *  gid = -1
 *  coes: zeros
 *  syms: zeros
 */
struct coded_packet *alloc_empty_packet(int size_g, int size_p)
{
	struct coded_packet *pkt = calloc(1, sizeof(struct coded_packet));
	if (pkt == NULL)
		return NULL;
	pkt->coes = calloc(size_g, sizeof(GF_ELEMENT));
	if (pkt->coes == NULL) 
		goto AllocErr;
	pkt->syms = calloc(size_p, sizeof(GF_ELEMENT));
	if (pkt->syms == NULL) 
		goto AllocErr;

	return pkt;

AllocErr:
	free_gnc_packet(pkt);
	return NULL;
}

/* Generate a GNC coded packet. Memory is allocated in the function. */
struct coded_packet *generate_gnc_packet(struct gnc_context *gc)
{
	struct coded_packet *pkt = alloc_empty_packet(gc->meta.size_g, gc->meta.size_p);
	int gid = schedule_generation(gc);
	encode_packet(gc, gid, pkt);
	return pkt;
}

/*
 * Generate a GNC coded packet in a given memory area.
 * It is the caller's responsibity to allocate memory properly.
 */
int generate_gnc_packet_im(struct gnc_context *gc, struct coded_packet *pkt)
{
	if (pkt == NULL || pkt->coes == NULL || pkt->syms == NULL)
		return -1;
	memset(pkt->coes, 0, gc->meta.size_g*sizeof(GF_ELEMENT));
	memset(pkt->syms, 0, gc->meta.size_p*sizeof(GF_ELEMENT));
	int gid = schedule_generation(gc);
	encode_packet(gc, gid, pkt);
	return (0);
}

void free_gnc_packet(struct coded_packet *pkt)
{
	if (pkt == NULL)
		return;
	if (pkt->coes != NULL)
		free(pkt->coes);
	if (pkt->syms != NULL)
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

/*
 * Print code summary
 * If called by decoders, it prints overhead and operations as well.
 */
void print_code_summary(struct gnc_metainfo *meta, int overhead, long long operations)
{
	char typestr[20];
	switch(meta->type) {
		case RAND_GNC_CODE:
			strcpy(typestr, "RAND");
			break;
		case BAND_GNC_CODE:
			strcpy(typestr, "BAND");
			break;
		default:
			strcpy(typestr, "UNKNOWN");
	}
	printf("datasize: %d ", meta->datasize);
	printf("precode: %.3f ", meta->pcrate);
	printf("size_b: %d ", meta->size_b);
	printf("size_g: %d ", meta->size_g);
	printf("size_p: %d ", meta->size_p);
	printf("type: %s ", typestr);
	printf("snum: %d ", meta->snum);
	printf("cnum: %d ", meta->cnum);
	printf("gnum: %d ", meta->gnum);
	if (operations != 0) {
		printf("overhead: %.3f ", (double) overhead/meta->snum);
		printf("computation: %f\n", (double) operations/meta->snum/meta->size_p);
	}
}

