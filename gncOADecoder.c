#include "gncOADecoder.h"

/*
 * Transform coefficient matrix of each generation to row echelon form (REF).
 * The REF has the property that if a diagonal element is nonzero, all elements
 * above it are reduced to zero.
 */
static long running_matrix_to_REF(struct decoding_context_OA *dec_ctx);

/*
 * Construct global decoding matrix (GDM) from REFs of each generation. Two times of
 * pivoting are done to transform GDM to a specific form where the bottom-right corner is
 * dense and the above half is partially diagonal:
 * -                 -
 * | x 0 0 0 0 x x x |
 * | 0 x 0 0 0 x x x |
 * | 0 0 x 0 0 x x x |
 * | 0 0 0 x 0 x x x |
 * | 0 0 0 0 x x x x |
 * | 0 0 0 0 0 x x x |
 * | 0 0 0 0 0 x x x |
 * | 0 0 0 0 0 x x x |
 * -                -
 */
static void construct_GDM(struct decoding_context_OA *dec_ctx);

/*
 * Fully transform GDM to identity matrix to finish decoding.
 */
static void diagonalize_GDM(struct decoding_context_OA *dec_ctx);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long pivot_matrix(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int *otoc, int *inactives);

/*
 * create_decoding_context_OA
 * Create context for overlap-aware (OA) decoding
 *  aoh - allowed overhead >=0
 */
void create_decoding_context_OA(struct decoding_context_OA *dec_ctx, long datasize, int s_b, int s_g, int s_p, int type, int aoh)
{
	static char fname[] = "create_decoding_context_OA";
	int i, j, k;

	// GNC code context
	// Since this is decoding, we construct GNC context without data
	// gc->pp will be filled by decoded packets
	struct gnc_context *gc;
	if (create_gnc_context(NULL, datasize, &gc, s_b, s_g, s_p, type) != 0) 
		printf("%s: create decoding context failed", fname);

	dec_ctx->gc = gc;

	dec_ctx->aoh		= aoh;
	dec_ctx->finished   = 0;
	dec_ctx->OA_ready	= 0;
	dec_ctx->local_DoF  = 0;
	dec_ctx->global_DoF = 0;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	// Allocate matrices for per-generation decoding
	dec_ctx->Matrices = calloc(dec_ctx->gc->meta.gnum, sizeof(struct running_matrix*));
	if (dec_ctx->Matrices == NULL)
		printf("%s: calloc dec_ctx->Matrices\n", fname);
	for (i=0; i<dec_ctx->gc->meta.gnum; i++) { 
		dec_ctx->Matrices[i] = calloc(1, sizeof(struct running_matrix));
		if (dec_ctx->Matrices[i] == NULL)
			printf("%s: malloc dec_ctx->Matrices[%d]\n", fname, i);
		dec_ctx->Matrices[i]->overhead = 0;
		
		// Allocate coefficient and message matrices in running_matrix
		// coefficeint: size_g x size_g
		// message:		size_g x size_p
		// Dim-1) Pointers to each row
		dec_ctx->Matrices[i]->coefficient = calloc(gensize, sizeof(GF_ELEMENT*));
		if (dec_ctx->Matrices[i]->coefficient == NULL)
			printf("%s: calloc dec_ctx->Matrices[%d]->coefficient\n", fname, i);
		dec_ctx->Matrices[i]->message = calloc(gensize, sizeof(GF_ELEMENT*));
		if (dec_ctx->Matrices[i]->message == NULL)
			printf("%s: calloc dec_ctx->Matrices[%d]->messsage\n", fname, i);
		// Dim-2) Elements of each row
		for (j=0; j<gensize; j++) {
			dec_ctx->Matrices[i]->coefficient[j] = calloc(gensize, sizeof(GF_ELEMENT));
			if (dec_ctx->Matrices[i]->coefficient[j] == NULL)
				printf("%s: calloc dec_ctx->Matrices[%d]->coefficient[%d]\n", fname, i, j);
			dec_ctx->Matrices[i]->message[j] = calloc(pktsize, sizeof(GF_ELEMENT));
			if (dec_ctx->Matrices[i]->message[j] == NULL)
				printf("%s: calloc dec_ctx->Matrices[%d]->messsage[%d]\n", fname, i, j);
		}
	}

	// Allocate matrices for global decoding
	dec_ctx->JMBcoefficient = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
	if (dec_ctx->JMBcoefficient == NULL)
		printf("%s: calloc dec_ctx->JMBcoefficient\n", fname);
	dec_ctx->JMBmessage     = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
	if (dec_ctx->JMBmessage == NULL)
		printf("%s: calloc dec_ctx->JMBmessage\n", fname);

	for (i=0; i<numpp+dec_ctx->aoh; i++) {
		dec_ctx->JMBcoefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
		dec_ctx->JMBmessage[i]     = calloc(pktsize, sizeof(GF_ELEMENT));
	}

	dec_ctx->inactives   = 0;
	dec_ctx->otoc_mapping = malloc(sizeof(int) * numpp);
	dec_ctx->ctoo_mapping = malloc(sizeof(int) * numpp);
	for (j=0; j<numpp; j++) {
		dec_ctx->otoc_mapping[j]   = j;				// original to current mapping
		dec_ctx->ctoo_mapping[j]   = j;				// current to original mapping
	}

	// JMB_coeffcient矩阵最下面的CHECKS行对应的是precoding generation，这部分从一开始就固定下来
	for (i=0; i<dec_ctx->gc->meta.cnum; i++) {
		dec_ctx->JMBcoefficient[dec_ctx->gc->meta.snum+dec_ctx->aoh+i][dec_ctx->gc->meta.snum+i] = 1;

		NBR_node *variable_node = dec_ctx->gc->graph->l_nbrs_of_r[i]->first; 		//ldpc_graph->nbrs_of_right[i];
		while (variable_node != NULL) {
			// 标记与该check packet连结的所有source packet node
			int src_pktid = variable_node->data;						//variable_node->nb_index;
			dec_ctx->JMBcoefficient[dec_ctx->gc->meta.snum+dec_ctx->aoh+i][src_pktid] = 1;
			//dec_ctx->JMBcoefficient[NUM_SRC+OHS+i][src_pktid] = variable_node->nb_ce;
			variable_node = variable_node->next;
		}
	}

	// performance indices
	dec_ctx->operations = 0;
	dec_ctx->overhead 	= 0;
}

void process_packet_OA(struct decoding_context_OA *dec_ctx, struct coded_packet *pkt)
{
	dec_ctx->overhead += 1;

	int i, j, k;
	int pivot;
	int pivotfound;
	GF_ELEMENT quotient;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	// Translate the encoding vector to the sorted form as in the generation
	int gid = pkt->gid;
	struct running_matrix *matrix = dec_ctx->Matrices[gid];
	matrix->overhead += 1;

	// start processing
	// Always process it in the generation first
	pivotfound = 0;
	for (i=0; i<gensize; i++) {
		if (pkt->coes[i] != 0) {
			if (matrix->coefficient[i][i] != 0) {
				quotient = galois_divide(pkt->coes[i], matrix->coefficient[i][i], GF_ORDER);
				dec_ctx->operations += 1;
				galois_multiply_add_region(&(pkt->coes[i]), &(matrix->coefficient[i][i]), quotient, gensize-i, GF_ORDER);
				dec_ctx->operations += (gensize - i);
				galois_multiply_add_region(pkt->syms, matrix->message[i], quotient, pktsize, GF_ORDER);
				dec_ctx->operations += pktsize;
			} else {
				pivotfound = 1;
				pivot = i;
				break;
			}
		}
	}
	
	/*
	 * If decoder is not OA ready, process the packet within the generation.
	 */
	if (dec_ctx->OA_ready != 1) {
		// cache as normal GNC packet
		if (pivotfound == 1) {
			memcpy(matrix->coefficient[pivot], pkt->coes, gensize*sizeof(GF_ELEMENT));
			memcpy(matrix->message[pivot], pkt->syms, pktsize*sizeof(GF_ELEMENT));
			dec_ctx->local_DoF += 1;
		}

		if ((dec_ctx->local_DoF >= dec_ctx->gc->meta.snum) && (dec_ctx->overhead >= (dec_ctx->gc->meta.snum+dec_ctx->aoh))) {
			dec_ctx->OA_ready = 1;
			// When OA ready, convert LDMs to upper triangular form
			long ops = running_matrix_to_REF(dec_ctx);
			dec_ctx->operations += ops;
			// Combine LDMs to GDM and apply inactivation pivoting
			construct_GDM(dec_ctx);

			// If numpp innovative packets are received, recover all
			// source packets from JMBcoeffcient and JMBmessage
			if (dec_ctx->global_DoF == numpp) 
				diagonalize_GDM(dec_ctx);
		}
	} else {
	/*
	 * If decoder is already OA ready, process the packet against the global matrix.
	 */
		if (pivotfound == 1) {
			// store into each generation anyway
			memcpy(matrix->coefficient[pivot], pkt->coes, gensize*sizeof(GF_ELEMENT));
			memcpy(matrix->message[pivot], pkt->syms,  pktsize*sizeof(GF_ELEMENT));

			/*
			 * Since the decoder is OA ready, translate the local encoded vector (LEV) 
			 * to global encoding vector (GEV). Since the GDM was probably pivoted, need
			 * to transform the GEV according to the pivoting order.
			 */
			GF_ELEMENT *re_ordered = calloc(numpp, sizeof(GF_ELEMENT));	
			for (i=0; i<gensize; i++) {
				/* obtain current index position of pktid */
				int curr_pos = dec_ctx->otoc_mapping[dec_ctx->gc->gene[gid]->pktid[i]];
				re_ordered[curr_pos] = pkt->coes[i];
			}
	
			/*
			 * Process the reordered GEV against GDM
			 */
			pivotfound = 0;								
			pivot = -1;
			for (int m=0; m<numpp; m++) {
				if (re_ordered[m] != 0) {
					if (dec_ctx->JMBcoefficient[m][m] != 0) {
						// mask the encoding vector and message over the JMB decoding matrix
						GF_ELEMENT quotient = galois_divide(re_ordered[m], dec_ctx->JMBcoefficient[m][m], GF_ORDER);
						dec_ctx->operations += 1;
						if ( m < (numpp - dec_ctx->inactives) ) {
							// Only needs to be multiply-and-add to the inactive part
							// this saves computation
							int maa_start = numpp - dec_ctx->inactives;		
							galois_multiply_add_region(re_ordered+maa_start, &(dec_ctx->JMBcoefficient[m][maa_start]), quotient, dec_ctx->inactives, GF_ORDER);
							dec_ctx->operations += dec_ctx->inactives;
							re_ordered[m] = 0;
						} else {
							galois_multiply_add_region(re_ordered+m, &(dec_ctx->JMBcoefficient[m][m]), quotient, numpp-m, GF_ORDER);
							dec_ctx->operations += (numpp - m);
						}
						galois_multiply_add_region(pkt->syms, dec_ctx->JMBmessage[m], quotient, pktsize, GF_ORDER);
						dec_ctx->operations += pktsize;
					} else {
						pivotfound = 1;
						pivot = m;
						break;
					}
				}
			}

			if (pivotfound == 1) {
				memcpy(dec_ctx->JMBcoefficient[pivot], re_ordered, numpp*sizeof(GF_ELEMENT));
				memcpy(dec_ctx->JMBmessage[pivot], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
				dec_ctx->global_DoF += 1;

				if (dec_ctx->global_DoF == numpp) {
					// recover all source from JMBcoeffcient & JMBmessage matrix
					diagonalize_GDM(dec_ctx);
				}
			}
			free(re_ordered);
		}
	}

	free_gnc_packet(pkt);
	pkt = NULL;
}

void free_decoding_context_OA(struct decoding_context_OA *dec_ctx)
{
	int i, j, k;
	for (i=0; i<dec_ctx->gc->meta.gnum; i++){
		// Free each decoding matrix
		for (j=0; j<dec_ctx->gc->meta.size_g; j++) {
			free(dec_ctx->Matrices[i]->coefficient[j]);
			free(dec_ctx->Matrices[i]->message[j]);
		}
		free(dec_ctx->Matrices[i]->coefficient);
		free(dec_ctx->Matrices[i]->message);
		free(dec_ctx->Matrices[i]);
	}
	free(dec_ctx->Matrices);
	for (j=0; j<dec_ctx->gc->meta.snum+dec_ctx->gc->meta.cnum+dec_ctx->aoh; j++) {
		free(dec_ctx->JMBcoefficient[j]);
		free(dec_ctx->JMBmessage[j]);
	}
	free(dec_ctx->JMBcoefficient);
	free(dec_ctx->JMBmessage);
	free(dec_ctx->otoc_mapping);
	free(dec_ctx->ctoo_mapping);
	free_gnc_context(dec_ctx->gc);
	free(dec_ctx);
	dec_ctx = NULL;
}

static void diagonalize_GDM(struct decoding_context_OA *dec_ctx)
{
	static char fname[] = "finish_recovering_inactivation";
	int i, j, k;
	int pos;
	int pktid;
	
	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	// Recover inactivated packets
	printf("Finishing decoding...\n");
	printf("Recovering \"inactive\" packets...\n");
	int ias = dec_ctx->inactives;
	GF_ELEMENT **ces_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
	for (i=0; i<ias; i++) {
		ces_submatrix[i] = calloc(ias, sizeof(GF_ELEMENT));
		memcpy(ces_submatrix[i], &(dec_ctx->JMBcoefficient[numpp-ias+i][numpp-ias]), ias*sizeof(GF_ELEMENT));
	}
	
	/* Perform back substitution to reduce the iasxias matrix to identity matrix */
	long long ops = back_substitute(ias, ias, pktsize, ces_submatrix, &(dec_ctx->JMBmessage[numpp-ias]));
	dec_ctx->operations += ops;

	// Recover decoded overlapping packets
	for (i=0; i<ias; i++) {
		// get original pktid at column (numpp-ias+i0
		pktid = dec_ctx->ctoo_mapping[numpp-ias+i];	
		// Construct decoded packets
		if ( (dec_ctx->gc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
			printf("%s: calloc gc->pp[%d]\n", fname, pktid);
		memcpy(dec_ctx->gc->pp[pktid], dec_ctx->JMBmessage[numpp-ias+i], sizeof(GF_ELEMENT)*pktsize);
	}
	/* free ces_submatrix */
	for (i=0; i<ias; i++) 
		free(ces_submatrix[i]);
	free(ces_submatrix);
	

	// Recover active packets
	printf("Recovering \"active\" packets...\n");
	GF_ELEMENT quotient;
	for (i=0; i<numpp-ias; i++) {
		/* 
		 * Clean up the inactive part of the upper half of GDM by
		 * masking non-zero element aginst already decoded inactive packets 
		 *
		 */
		for (j=numpp-ias; j<numpp; j++) {
			if (dec_ctx->JMBcoefficient[i][j] != 0) {
				quotient = dec_ctx->JMBcoefficient[i][j];
				pktid = dec_ctx->ctoo_mapping[j];			
				galois_multiply_add_region(dec_ctx->JMBmessage[i], dec_ctx->gc->pp[pktid], quotient, pktsize, GF_ORDER);
				dec_ctx->JMBcoefficient[i][j] = 0;
				dec_ctx->operations += pktsize;
			}
		}

		// Convert diagonal elements of top-left part of T to 1
		quotient = dec_ctx->JMBcoefficient[i][i];
		if (quotient != 1) {
			for (int n=0; n<pktsize; n++)
				dec_ctx->JMBmessage[i][n] = galois_divide(dec_ctx->JMBmessage[i][n], quotient, GF_ORDER);
			dec_ctx->operations += pktsize;
			dec_ctx->JMBcoefficient[i][i] = 1;
		}

		// Save the decoded packet
		pktid = dec_ctx->ctoo_mapping[i];
		if ( dec_ctx->gc->pp[pktid] != NULL )
			printf("%s：warning: packet %d is already recovered.\n", fname, pktid);
		if ( (dec_ctx->gc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
			printf("%s: calloc gc->pp[%d]\n", fname, pktid);
		memcpy(dec_ctx->gc->pp[pktid], dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT)*pktsize);
	}

	dec_ctx->finished = 1;
}

// Partially diagonalize all running matrices when the decoder is OA ready
// diagonalizes them as much as possible
static long running_matrix_to_REF(struct decoding_context_OA *dec_ctx)
{
	long long operations = 0;
	int i, j, k, l;
	GF_ELEMENT quotient;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;


	for (i=0; i<dec_ctx->gc->meta.gnum; i++) {
		struct running_matrix *matrix = dec_ctx->Matrices[i];

		// Partially diagonalize the LDM
		int consecutive = 1;
		for (k=gensize-1; k>=0; k--) {
			if (matrix->coefficient[k][k] == 0) {
				consecutive = 0;
				continue;
			}

			// eliminate elements above the nonzero diagonal elements
			for (l=0; l<k; l++) {
				if (matrix->coefficient[l][k] == 0)
					continue;
				
				quotient = galois_divide(matrix->coefficient[l][k], matrix->coefficient[k][k], GF_ORDER);
				operations += 1;
				matrix->coefficient[l][k] = 0;
				// 注意后面有的列可能并非为零列(也就是对角元素非零)，这些也要eliminate
				for (int m=k+1; m<gensize; m++) {
					if (matrix->coefficient[m][m] == 0) {
						matrix->coefficient[l][m] = galois_add(matrix->coefficient[l][m], galois_multiply(matrix->coefficient[k][m], quotient, GF_ORDER));
						operations += 1;
					}
				}
				galois_multiply_add_region(matrix->message[l], matrix->message[k], quotient, pktsize, GF_ORDER);
				operations += pktsize;
			}
		}
		if (consecutive == 1) 
			printf("Class %d is self-decodable.\n", i);
	}
	return operations;
}

static void construct_GDM(struct decoding_context_OA *dec_ctx)
{
	int i, j, k;
	struct running_matrix *matrix;

	int gensize = dec_ctx->gc->meta.size_g;
	int pktsize = dec_ctx->gc->meta.size_p;
	int numpp   = dec_ctx->gc->meta.snum + dec_ctx->gc->meta.cnum;

	// Step 1, translate LEVs to GEV and move them to GDM
	GF_ELEMENT *global_ces = calloc(numpp, sizeof(GF_ELEMENT));
	int p_copy = 0;								// 拷贝到JMBcofficient的行指针
	for (i=0; i<dec_ctx->gc->meta.gnum; i++) {
		matrix = dec_ctx->Matrices[i];
		for (j=0; j<gensize; j++) {
			if (matrix->coefficient[j][j] == 0)
				continue;						// there is no local DoF here
			else {
				memset(global_ces, 0, numpp*sizeof(GF_ELEMENT));	/* Reset before reuse */
				for (k=0; k<gensize; k++) {
					//global_ces[matrix->indices[k]] = matrix->coefficient[j][k];
					global_ces[dec_ctx->gc->gene[i]->pktid[k]] = matrix->coefficient[j][k];
				}
				memcpy(dec_ctx->JMBcoefficient[p_copy], global_ces, numpp*sizeof(GF_ELEMENT));
				memcpy(dec_ctx->JMBmessage[p_copy], matrix->message[j], pktsize*sizeof(GF_ELEMENT));
				p_copy += 1;
			}
		}
	}
	free(global_ces);
	printf("%d local DoFs are available, copied %d to GDM.\n", dec_ctx->local_DoF, p_copy);

	/* Transform GDM using 2 rounds of pivoting */
	//dec_ctx->operations += pivot_matrix(numpp+dec_ctx->aoh, numpp, pktsize, dec_ctx->JMBcoefficient, dec_ctx->JMBmessage, dec_ctx->otoc_mapping, dec_ctx->ctoo_mapping, &(dec_ctx->inactives));
	dec_ctx->operations += pivot_matrix(numpp+dec_ctx->aoh, numpp, pktsize, dec_ctx->JMBcoefficient, dec_ctx->JMBmessage, dec_ctx->otoc_mapping, &(dec_ctx->inactives));
	/* Count available innovative rows */
	for (i=0; i<numpp; i++) {
		if (dec_ctx->JMBcoefficient[i][i] != 0)
			dec_ctx->global_DoF++;
		dec_ctx->ctoo_mapping[dec_ctx->otoc_mapping[i]] = i;
	}
}
