/********************************************************************
 *                       Compact Band Decoder
 *
 * The decoder only applies to band GNC code. 
 *
 * Unlike regular band decoder (gncBandDecoder.c), CBD decoder is 
 * compact in coefficient matrix storage. Only coefficients in the band
 * (which are nonzeros) are stored. A price to pay is that pivoting 
 * cannot be performed due to the limited random access and row/col 
 * manipulation capability of using compact row vectors.
 ********************************************************************/
#include "common.h"
#include "galois.h"
#include "bipartite.h"
#include "gncCBDDecoder.h"
static int process_vector_CBD(struct decoding_context_CBD *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message);
static int apply_parity_check_matrix(struct decoding_context_CBD *dec_ctx);
static void finish_recovering_CBD(struct decoding_context_CBD *dec_ctx);

// create decoding context for band decoder
void create_decoding_context_CBD(struct decoding_context_CBD *dec_ctx, long datasize, struct gnc_parameter gp)
{
	static char fname[] = "create_decoding_context_CBD";
	int i, j, k;

	// GNC code context
	// Since this is decoding, we construct GNC context without data
	// gc->pp will be filled by decoded packets
	if (gp.type != BAND_GNC_CODE) {
		fprintf(stderr, "Band decoder only applies to band GNC code.\n");
		return;
	}
	struct gnc_context *gc;
	if (create_gnc_context(NULL, datasize, &gc, gp) != 0) 
		fprintf(stderr, "%s: create decoding context failed", fname);

	dec_ctx->gc = gc;

	dec_ctx->finished     = 0;
	dec_ctx->DoF 		  = 0;
	dec_ctx->de_precode   = 0;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	dec_ctx->row = (struct row_vector **) calloc(numpp, sizeof(struct row_vector *));
	if (dec_ctx->row == NULL)
		fprintf(stderr, "%s: calloc dec_ctx->row failed\n", fname);
	dec_ctx->message = calloc(numpp, sizeof(GF_ELEMENT*));
	if (dec_ctx->message == NULL)
		fprintf(stderr, "%s: calloc dec_ctx->message failed\n", fname);
	for (i=0; i<numpp; i++) { 
		dec_ctx->message[i] = calloc(pktsize, sizeof(GF_ELEMENT));
		if (dec_ctx->message[i] == NULL)
			fprintf(stderr, "%s: calloc dec_ctx->message[%d] failed\n", fname, i);
	}

	dec_ctx->overhead     = 0;
	dec_ctx->operations   = 0;
}

/* 
 * Note: throughout the packet collecting process, the decoding matrix 
 * is maintained an upper triangular form.
 */
void process_packet_CBD(struct decoding_context_CBD *dec_ctx, struct coded_packet *pkt)
{
	static char fname[] = "process_packet_CBD";
	dec_ctx->overhead += 1;
	int i, j, k;
	GF_ELEMENT quotient;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	// transform GNC encoding vector to full length
	GF_ELEMENT *ces = calloc(numpp, sizeof(GF_ELEMENT));
	if (ces == NULL)
		fprintf(stderr, "%s: calloc ces failed\n", fname);
	for (i=0; i<gensize; i++) {
		int index = dec_ctx->gc->gene[pkt->gid]->pktid[i];
		ces[index] = pkt->coes[i];
	}

	/* Process full-length encoding vector against decoding matrix */
	int pivot = process_vector_CBD(dec_ctx, ces, pkt->syms);
	free(ces);
	ces = NULL;
	free_gnc_packet(pkt);
	pkt = NULL;
	// If the number of received DoF is equal to NUM_SRC, apply the parity-check matrix.
	// The messages corresponding to rows of parity-check matrix are all-zero.
	
	if (dec_ctx->DoF == dec_ctx->gc->meta.snum) {
		dec_ctx->de_precode = 1;	/*Mark de_precode before applying precode matrix*/
		int missing_DoF = apply_parity_check_matrix(dec_ctx);
#if defined(GNCTRACE)
		printf("After applying the parity-check matrix, %d DoF are missing.\n", missing_DoF);
#endif
		dec_ctx->DoF = numpp - missing_DoF;
	}

	if (dec_ctx->DoF == dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum) {
		finish_recovering_CBD(dec_ctx);
	}
}

/* Process a full row vector against CBD decoding matrix */
static int process_vector_CBD(struct decoding_context_CBD *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message)
{
	static char fname[] = "process_vector_CBD";
	int i, j, k;
	int pivot = -1;
	int pivotfound = 0;
	GF_ELEMENT quotient;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	for (i=0; i<numpp; i++) {
		if (vector[i] != 0) {
			if (dec_ctx->row[i] != NULL) {
			/* There is a valid row saved for pivot-i, process against it */
				assert(dec_ctx->row[i]->elem[0]);
				quotient = galois_divide(vector[i], dec_ctx->row[i]->elem[0], GF_POWER);
				dec_ctx->operations += 1;
				galois_multiply_add_region(&(vector[i]), dec_ctx->row[i]->elem, quotient, dec_ctx->row[i]->len, GF_POWER);
				assert(!vector[i]);
				dec_ctx->operations += dec_ctx->row[i]->len;
				galois_multiply_add_region(message, dec_ctx->message[i], quotient, pktsize, GF_POWER);
				dec_ctx->operations += pktsize;
			} else {
				pivotfound = 1;
				pivot = i;
				break;
			}
		}
	}

	if (pivotfound == 1) {
		/* Save it to the corresponding row */
		dec_ctx->row[pivot] = (struct row_vector*) malloc(sizeof(struct row_vector));
		if (dec_ctx->row[pivot] == NULL)
			fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, pivot);
		int len;
		if (!dec_ctx->de_precode) { 
			/* before de_precode every row is no more than gensize-width */
			len = numpp - pivot > gensize ? gensize : numpp - pivot;
		} else {
			/* row bandwidth is indetermined, so being conservative here */
			len = numpp - pivot;
		}
		dec_ctx->row[pivot]->len = len;
		dec_ctx->row[pivot]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
		if (dec_ctx->row[pivot]->elem == NULL)
			fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, pivot);
		memcpy(dec_ctx->row[pivot]->elem, &(vector[pivot]), len*sizeof(GF_ELEMENT));
		assert(dec_ctx->row[pivot]->elem[0]);
		memcpy(dec_ctx->message[pivot], message,  pktsize*sizeof(GF_ELEMENT));
		dec_ctx->DoF += 1;
	}
	return pivot;
}

// Apply the parity-check matrix to the decoding matrix
static int apply_parity_check_matrix(struct decoding_context_CBD *dec_ctx)
{
	static char fname[] = "apply_parity_check_matrix";
#if defined(GNCTRACE)
	printf("Entering %s\n", fname);
#endif
	int i, j, k;
	int num_of_new_DoF = 0;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	// 1, Copy parity-check vectors to the nonzero rows of the decoding matrix
	GF_ELEMENT *ces = malloc(numpp*sizeof(GF_ELEMENT));
	GF_ELEMENT *msg = malloc(pktsize*sizeof(GF_ELEMENT));
	int p = 0;			// index pointer to the parity-check vector that is to be copyed
	for (int p=0; p<dec_ctx->gc->meta.cnum; p++) {
		memset(ces, 0, numpp*sizeof(GF_ELEMENT));
		memset(msg, 0, pktsize*sizeof(GF_ELEMENT));
		/* Set the coding vector according to parity-check bits */
		NBR_node *varnode = dec_ctx->gc->graph->l_nbrs_of_r[p]->first;
		while (varnode != NULL) {
			ces[varnode->data] = 1;
			varnode = varnode->next;
		}
		ces[dec_ctx->gc->meta.snum+p] = 1;
		int pivot = process_vector_CBD(dec_ctx, ces, msg);
	}
	free(ces);
	free(msg);

	/* Count available innovative rows */
	int missing_DoF = 0;
	for (i=0; i<numpp; i++) {
		if (dec_ctx->row[i] == NULL) 
			missing_DoF++;
#if defined(GNCTRACE)
		else if (dec_ctx->row[i]->elem[0] ==0)
			printf("%s: row[%d]->elem[0] is 0\n", fname, i);
#endif
	}
	return missing_DoF;
}


static void finish_recovering_CBD(struct decoding_context_CBD *dec_ctx)
{
	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;
	int i, j;
	int len;
	GF_ELEMENT quotient;
	for (i=numpp-1; i>=0; i--) {
		/* eliminate all nonzeros above diagonal elements from right to left*/
		for (j=0; j<i; j++) {
			len = dec_ctx->row[j]->len;
			if (j+len <= i || dec_ctx->row[j]->elem[i-j] == 0)
				continue;
			assert(dec_ctx->row[i]->elem[0]);
			quotient = galois_divide(dec_ctx->row[j]->elem[i-j], dec_ctx->row[i]->elem[0], GF_POWER);
			galois_multiply_add_region(dec_ctx->message[j], dec_ctx->message[i], quotient, pktsize, GF_POWER);
			dec_ctx->operations += (pktsize + 1);
			dec_ctx->row[j]->elem[i-j] = 0;
		}
		/* convert diagonal to 1*/
		if (dec_ctx->row[i]->elem[0] != 1) {
			galois_multiply_region(dec_ctx->message[i], galois_divide(1, dec_ctx->row[i]->elem[0], GF_POWER), pktsize, GF_POWER);
			dec_ctx->operations += (pktsize + 1);
			dec_ctx->row[i]->elem[0] = 1;
		}
		/* save decoded packet */
		dec_ctx->gc->pp[i] = calloc(pktsize, sizeof(GF_ELEMENT));
		memcpy(dec_ctx->gc->pp[i], dec_ctx->message[i], pktsize*sizeof(GF_ELEMENT));
	}
	dec_ctx->finished = 1;
}

void free_decoding_context_CBD(struct decoding_context_CBD *dec_ctx)
{
	for (int i=dec_ctx->gc->meta.snum+dec_ctx->gc->meta.cnum-1; i>=0; i--) {
		free(dec_ctx->row[i]->elem);
		free(dec_ctx->row[i]);
		free(dec_ctx->message[i]);
	}
	free(dec_ctx->row);
	free(dec_ctx->message);
	free_gnc_context(dec_ctx->gc);
	free(dec_ctx);
	dec_ctx = NULL;
}

