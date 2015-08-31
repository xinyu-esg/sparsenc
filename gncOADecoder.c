#include "gncOADecoder.h"

#define ZLATEVS			3				// the number of searching rows using Zlatev's strategy 
#define IA_INIT			0
#define IA_STEP			1				// gradually inactivate more columns
typedef struct subscript				Subscript;
typedef struct subscripts				ssList;

static void finish_recovering_inactivation(struct decoding_context_OA *dec_ctx);
static long long partially_diag_RM_matrices(struct decoding_context_OA *dec_ctx);
static void construct_GDM_inactivation(struct decoding_context_OA *dec_ctx);

static int Inactivation_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots);
static int Zlatev_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots);

/*
 * Reorder matrix A and B of Ax=B according to pivot sequence in (RowPivots, ColPivots)
 */
static void matrices_reordering(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, ssList *RowPivots, ssList *ColPivots, int npv);

/*
 * Reorder columns of A according to sequence in ColPivots
 */
static void permute_matrix_columns(int nrow, int ncolA, GF_ELEMENT **A, ssList *ColPivots);

/* 
 * Helper functions in handling pivots
 * We use double-linked lists when pivoting a matrix.
 */
static void insertSubAtBeginning(ssList *sub_list, Subscript *sub);
static void insertSubAtEnd(ssList *sub_list, Subscript *sub);
static void removeSubscript(ssList *sub_list, Subscript *sub);
static void free_subscriptList(ssList *sub_list);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);

/*
 * create_decoding_context_OA
 * Create context for overlap-aware (OA) decoding
 *  aoh - allowed overhead >=0
 */
void create_decoding_context_OA(struct decoding_context_OA *dec_ctx, long datasize, int s_b, int s_g, int s_p, int aoh)
{
	static char fname[] = "create_decoding_context_OA";
	int i, j, k;

	// GNC code context
	// Since this is decoding, we construct GNC context without data
	// gc->pp will be filled by decoded packets
	struct gnc_context *gc;
	if (create_gnc_context(NULL, datasize, &gc, s_b, s_g, s_p) != 0) 
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
			int src_pkt_id = variable_node->data;						//variable_node->nb_index;
			dec_ctx->JMBcoefficient[dec_ctx->gc->meta.snum+dec_ctx->aoh+i][src_pkt_id] = 1;
			//dec_ctx->JMBcoefficient[NUM_SRC+OHS+i][src_pkt_id] = variable_node->nb_ce;
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
	GF_ELEMENT *local_ces = calloc(gensize, sizeof(GF_ELEMENT));		//[matrix->degree];
	memcpy(local_ces, pkt->coes, sizeof(GF_ELEMENT)*gensize);

	// start processing
	// Always process it in the generation first
	pivotfound = 0;
	for (i=0; i<gensize; i++) {
		if (local_ces[i] != 0) {
			if (matrix->coefficient[i][i] != 0) {
				quotient = galois_divide(local_ces[i], matrix->coefficient[i][i], GF_ORDER);
				dec_ctx->operations += 1;
				galois_multiply_add_region(local_ces+i, &(matrix->coefficient[i][i]), quotient, gensize-i, GF_ORDER);
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
			memcpy(matrix->coefficient[pivot], local_ces, gensize*sizeof(GF_ELEMENT));
			memcpy(matrix->message[pivot], pkt->syms, pktsize*sizeof(GF_ELEMENT));
			dec_ctx->local_DoF += 1;
		}

		if ((dec_ctx->local_DoF >= dec_ctx->gc->meta.snum) && (dec_ctx->overhead >= (dec_ctx->gc->meta.snum+dec_ctx->aoh))) {
			dec_ctx->OA_ready = 1;
			// When OA ready, convert LDMs to upper triangular form
			long long ops = partially_diag_RM_matrices(dec_ctx);
			dec_ctx->operations += ops;
			// Combine LDMs to GDM and apply inactivation pivoting
			construct_GDM_inactivation(dec_ctx);

			// If numpp innovative packets are received, recover all
			// source packets from JMBcoeffcient and JMBmessage
			if (dec_ctx->global_DoF == numpp) 
				finish_recovering_inactivation(dec_ctx);
		}
	} else {
	/*
	 * If decoder is already OA ready, process the packet against the global matrix.
	 */
		if (pivotfound == 1) {
			// store into each generation anyway
			memcpy(matrix->coefficient[pivot], local_ces, gensize*sizeof(GF_ELEMENT));
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
				re_ordered[curr_pos] = local_ces[i];
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
					finish_recovering_inactivation(dec_ctx);
				}
			}
			free(re_ordered);
		}
	}

	free(local_ces);			// NOTE: is it really necessary to allocate local_ces?
	free_gnc_packet(pkt);
	pkt = NULL;
}

static void finish_recovering_inactivation(struct decoding_context_OA *dec_ctx)
{
	static char fname[] = "finish_recovering_inactivation";
	int i, j, k;
	int pos;
	int pkt_id;
	
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
		// NOTE: here better to use memcpy
		for (j=0; j<ias; j++) 
			ces_submatrix[i][j] = dec_ctx->JMBcoefficient[numpp-ias+i][numpp-ias+j];
	}
	GF_ELEMENT **msg_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
	for (i=0; i<ias; i++) {
		msg_submatrix[i] = calloc(pktsize, sizeof(GF_ELEMENT));
		// NOTE: here better to use memcpy
		memcpy(msg_submatrix[i], dec_ctx->JMBmessage[numpp-ias+i], sizeof(GF_ELEMENT)*pktsize);
	}

	//long long ops = back_substitute(ias, ias, dec_ctx->gc->meta.size_p, ces_submatrix, dec_ctx->JMBmessage+dec_ctx->gc->meta.snum+dec_ctx->gc->meta.cnum-ias);
	//long long ops = back_substitute(ias, ias, pktsize, ces_submatrix, &(dec_ctx->JMBmessage[numpp-ias]));
	long long ops = back_substitute(ias, ias, pktsize, ces_submatrix, msg_submatrix);
	// copy msg_submatrix back to JMBmessage
	for (i=0; i<ias; i++) 
		memcpy(dec_ctx->JMBmessage[numpp-ias+i], msg_submatrix[i], sizeof(GF_ELEMENT)*pktsize);
	dec_ctx->operations += ops;

	// Recover decoded overlapping packets
	for (i=0; i<ias; i++) {
		// get original pkt_id at column (numpp-ias+i0
		int pkt_id = dec_ctx->ctoo_mapping[numpp-ias+i];	
		// Construct decoded packets
		if (dec_ctx->gc->pp[pkt_id] != NULL)
			printf("%s: Warning: dec_ctx->gc->pp[%d] is already there.\n", fname, pkt_id);
		if ( (dec_ctx->gc->pp[pkt_id] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
			printf("%s: calloc gc->pp[%d]\n", fname, pkt_id);
		memcpy(dec_ctx->gc->pp[pkt_id], dec_ctx->JMBmessage[numpp-ias+i], sizeof(GF_ELEMENT)*pktsize);
	}
	// free ces_submatrix
	for (i=0; i<ias; i++) 
		free(ces_submatrix[i]);
	free(ces_submatrix);
	// free msg_submatrix
	for (i=0; i<ias; i++) 
		free(msg_submatrix[i]);
	free(msg_submatrix);
	

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
				int pkt_id = dec_ctx->ctoo_mapping[j];			
				if (dec_ctx->gc->pp[pkt_id] == NULL)
					printf("%s: error: packet %d is not decoded yet.\n", fname, pkt_id);
				galois_multiply_add_region(dec_ctx->JMBmessage[i], dec_ctx->gc->pp[pkt_id], quotient, pktsize, GF_ORDER);
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
		int pkt_id = dec_ctx->ctoo_mapping[i];
		if ( dec_ctx->gc->pp[pkt_id] != NULL )
			printf("%s：warning: packet %d is already recovered.\n", fname, pkt_id);
		if ( (dec_ctx->gc->pp[pkt_id] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
			printf("%s: calloc gc->pp[%d]\n", fname, pkt_id);
		memcpy(dec_ctx->gc->pp[pkt_id], dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT)*pktsize);
	}

	dec_ctx->finished = 1;
}

// Partially diagonalize all running matrices when the decoder is OA ready
// diagonalizes them as much as possible
static long long partially_diag_RM_matrices(struct decoding_context_OA *dec_ctx)
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

static void construct_GDM_inactivation(struct decoding_context_OA *dec_ctx)
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

	// Step 2, pivoting the GDM using inactivation pivoting
	ssList *row_pivots = malloc(sizeof(ssList));
	ssList *col_pivots = malloc(sizeof(ssList));
	row_pivots->ssFirst = row_pivots->ssLast = NULL;
	col_pivots->ssFirst = col_pivots->ssLast = NULL;
	double pivoting_time=0.0;
	time_t start_pivoting, stop_pivoting;
	printf("Start pivoting...\n");
	time(&start_pivoting);
	
	// 拷贝一份GDM和message matrix用于pivoting及后面的处理
	GF_ELEMENT **ces_matrix = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
	GF_ELEMENT **msg_matrix = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
	for (i=0; i<numpp+dec_ctx->aoh; i++) {
		ces_matrix[i] = calloc(numpp, sizeof(GF_ELEMENT));
		memcpy(ces_matrix[i], dec_ctx->JMBcoefficient[i], numpp * sizeof(GF_ELEMENT));
		msg_matrix[i] = calloc(pktsize, sizeof(GF_ELEMENT));
		memcpy(msg_matrix[i], dec_ctx->JMBmessage[i], pktsize * sizeof(GF_ELEMENT));
	}

	int npv;				// the number of pivots found (这是个无用的量,placeholder)
	int ias;				// the number of columns been inactivated
	printf("Inactivation pivoting strategy is used. IA_INIT: %d, IA_STEP: %d.\n", IA_INIT, IA_STEP);
	ias = Inactivation_pivoting(numpp+dec_ctx->aoh, numpp, ces_matrix, row_pivots, col_pivots);
	dec_ctx->inactives = ias;
	printf("A total of %d/%d columns are inactivated.\n", ias, numpp);
	npv = numpp;
	time(&stop_pivoting);
	pivoting_time = difftime(stop_pivoting, start_pivoting);
	printf("Time consumed in pivoting T: %.0f seconds.\n", pivoting_time);

	// 保存re-order过后列的顺序
	// Inactivation后第一次re-ordering，记下indices的mapping
	Subscript *ss_pt = col_pivots->ssFirst;
	for (i=0; i<npv; i++) {
		dec_ctx->otoc_mapping[ss_pt->index] = i;		// 第i个pivot对应的原来第ss_pt->index列
		dec_ctx->ctoo_mapping[i] = ss_pt->index;	
		ss_pt = ss_pt->next;
	}

	printf("Start matrix re-ordering after the inactivation...\n");
	double reordering_time = 0.0;
	time_t start_reorder, stop_reorder;
	time(&start_reorder);
	matrices_reordering(numpp+dec_ctx->aoh, numpp, pktsize, ces_matrix, msg_matrix, row_pivots, col_pivots, npv);
	time(&stop_reorder);
	reordering_time = difftime(stop_reorder, start_reorder);
	printf("Time consumed in matrix re-ordering after pivoting T: %.0f seconds.\n", reordering_time);

	free_subscriptList(row_pivots);
	free_subscriptList(col_pivots);

	int missing_pivots = 0;
	for (i=0; i<numpp+dec_ctx->aoh; i++) {
		if (ces_matrix[i][i] == 0) {
			missing_pivots += 1;
		}
	}
	printf("There are %d pivots missing\n", missing_pivots);

	/********
	 *
	 * Check the density of the lower half of the active part
	 *
	 *******/
	
	long long nonzeros = 0;
	for (i=numpp-ias; i<numpp; i++) {
		for (j=0; j<numpp-ias; j++)
			if (ces_matrix[i][j] != 0)
				nonzeros++;
	}
	printf("Fraction of nonzeros in the lower half of the active part: %f\n", (double) nonzeros/(numpp-ias)/ias);


	//将active part变为对角阵
	printf("Convert the left half of T (lower triangular) to diagonal...\n");	
	long long ops1=0;
	GF_ELEMENT quotient;
	for (i=0; i<numpp-ias; i++) {
		for (j=i+1; j<numpp+dec_ctx->aoh; j++) {
			if (ces_matrix[i][i] == 0)
				printf("The diagonal element after re-ordering is nonzero.\n");

			// process the item on (j, i)
			if (ces_matrix[j][i] != 0) {
				quotient = galois_divide(ces_matrix[j][i], ces_matrix[i][i], GF_ORDER);
				ops1 += 1;
				// XOR the corresponding part in the inactive part
				galois_multiply_add_region(&(ces_matrix[j][numpp-ias]), &(ces_matrix[i][numpp-ias]), quotient, ias, GF_ORDER);
				ops1 += ias;	// This part of matrix in processing is lower triangular in part, so operations only needed in the back
				// simultaneously do the same thing on right matrix B
				galois_multiply_add_region(msg_matrix[j], msg_matrix[i], quotient, pktsize, GF_ORDER);
				ops1 += pktsize;				
				ces_matrix[j][i] = 0;			// eliminate the item
			}
		}
	}
	dec_ctx->operations += ops1;
	
	/********
	 *
	 * Check the density of T_I
	 *
	 *******/
	/*
	 *
	int row_counts_tmp;
	for (i=nrow-ncol+ias; i<nrow; i++)
	{
		row_counts_tmp = 0;
		for (j=ncol-ias; j<ncol; j++)
		{
			if (ces_matrix[i][j] != 0)
				row_counts_tmp += 1;
		}
		printf("row %d has %d/%d nonzero elements.\n", i, row_counts_tmp, ias);
	}
	*/

	GF_ELEMENT **ces_submatrix = calloc(ias+dec_ctx->aoh, sizeof(GF_ELEMENT*));
	GF_ELEMENT **msg_submatrix = calloc(ias+dec_ctx->aoh, sizeof(GF_ELEMENT*));
	for (i=0; i<ias+dec_ctx->aoh; i++){
		ces_submatrix[i] = calloc(ias, sizeof(GF_ELEMENT));
		memcpy(ces_submatrix[i], &(ces_matrix[numpp-ias+i][numpp-ias]), ias*sizeof(GF_ELEMENT));
		msg_submatrix[i] = calloc(pktsize, sizeof(GF_ELEMENT));
		memcpy(msg_submatrix[i], msg_matrix[numpp-ias+i], pktsize*sizeof(GF_ELEMENT));
	}

	// 对inactive part的右下角的(ias x ias)的方阵再做一次pivoting
	printf("Perform another time of pivoting after partly diagonalizing T\n");
	ssList *row_pivots_2nd = malloc(sizeof(ssList));
	ssList *col_pivots_2nd = malloc(sizeof(ssList));
	row_pivots_2nd->ssFirst = row_pivots_2nd->ssLast = NULL;
	col_pivots_2nd->ssFirst = col_pivots_2nd->ssLast = NULL;
	printf("Start the second-time pivoting...\n");
	int ias_2nd = Zlatev_pivoting(ias+dec_ctx->aoh, ias, ces_submatrix, row_pivots_2nd, col_pivots_2nd);
	printf("A total of %d/%d columns are further inactivated.\n", ias_2nd, ias);
	int npv_2nd = ias;
	
	// 更新otoc_mapping & ctoo_mapping after second-time pivoting
	// Careful!!!这一步非常重要，也非常容易出错！
	int *partial_mappings = (int *) calloc(ias, sizeof(int));

	Subscript *ss_pt_2nd = col_pivots_2nd->ssFirst;
	for (i=0; i<npv_2nd; i++) {
		partial_mappings[i] = dec_ctx->ctoo_mapping[numpp-ias+ss_pt_2nd->index];
		ss_pt_2nd = ss_pt_2nd->next;
	}
	memcpy(&(dec_ctx->ctoo_mapping[numpp-ias]), partial_mappings, sizeof(int)*ias);
	free(partial_mappings);
	for (i=0; i<numpp; i++) {
		dec_ctx->otoc_mapping[dec_ctx->ctoo_mapping[i]] = i;
	}

	// re-order the ((ias+OHS) x ias) matrix
	matrices_reordering(ias+dec_ctx->aoh, ias, pktsize, ces_submatrix, msg_submatrix, row_pivots_2nd, col_pivots_2nd, npv_2nd);
	
	// re-order the columns above it accordingly
	GF_ELEMENT **BI_matrix = calloc(numpp+dec_ctx->aoh-ias, sizeof(GF_ELEMENT*));
	for (i=0; i<numpp+dec_ctx->aoh-ias; i++) {
		BI_matrix[i] = calloc(ias, sizeof(GF_ELEMENT));
		for (j=0; j<ias; j++) {
			BI_matrix[i][j] = ces_matrix[i][numpp-ias+j];
		}
	}

	permute_matrix_columns(numpp+dec_ctx->aoh-ias, ias, BI_matrix, col_pivots_2nd);
	// copy back
	for (i=0; i<numpp+dec_ctx->aoh-ias; i++) {
		for (j=0; j<ias; j++) {
			ces_matrix[i][numpp-ias+j] = BI_matrix[i][j];
		}
	}
	free_subscriptList(row_pivots_2nd);
	free_subscriptList(col_pivots_2nd);
	// 结束第二次pivoting

	// 3, perform forward substitution on the re-ordered matrix
	// Feb01注：此时只需要对T的右下角sub-matrix作Gaussian Elimination
	// 动态分配decoding matrix和message matrix，用来保存T右下角的矩阵及其对应的message矩阵
	printf("Perform forward substitution on the bottom-right part of GDM (size: %d x %d)...\n", ias, ias);
	long long ops = forward_substitute(ias+dec_ctx->aoh, ias, pktsize, ces_submatrix, msg_submatrix);
	dec_ctx->operations += ops;

	int missing_DoF = 0;
	for (i=0; i<ias; i++) {
		if (ces_submatrix[i][i] == 0) {
			missing_DoF += 1;
		}
	}
	printf("There are %d DoFs missing after performed forward substitution.\n", missing_DoF);
	dec_ctx->global_DoF = (numpp - missing_DoF);

	// store the processed matrix back to dec_ctx
	for (i=0; i<numpp-ias; i++) {
		for (j=0; j<numpp; j++) {
			dec_ctx->JMBcoefficient[i][j] = ces_matrix[i][j];
		}
		memcpy(dec_ctx->JMBmessage[i], msg_matrix[i], pktsize*sizeof(GF_ELEMENT));
	}

	for (i=numpp-ias; i<numpp; i++) {
		memset(dec_ctx->JMBcoefficient[i], 0, sizeof(GF_ELEMENT)*numpp);
		for (j=numpp-ias; j<numpp; j++) {
			dec_ctx->JMBcoefficient[i][j] = ces_submatrix[i-(numpp-ias)][j-(numpp-ias)];
		}
		memcpy(dec_ctx->JMBmessage[i], msg_submatrix[i-(numpp-ias)], pktsize*sizeof(GF_ELEMENT));
	}

	// free ces_matrix and msg_matrix
	for (i=0; i<numpp+dec_ctx->aoh; i++) {
		free(ces_matrix[i]);
		free(msg_matrix[i]);
	}
	free(ces_matrix);
	free(msg_matrix);

	//free ces_submatrix and msg_submatrix
	for (i=0; i<ias+dec_ctx->aoh; i++){
		free(ces_submatrix[i]);
		free(msg_submatrix[i]);
	}
	free(ces_submatrix);
	free(msg_submatrix);

	// free BI_matrix
	for (i=0; i<numpp+dec_ctx->aoh-ias; i++) {
		free(BI_matrix[i]);
	}
	free(BI_matrix);
}

// Use Zlatev pivoting scheme
static int Zlatev_pivoting(int nrow, int ncolA, GF_ELEMENT *A[], ssList *RowPivots, ssList *ColPivots)
{	
	// 对A的各行各列的非零元素计数
	//对A的各行各列的非零元素计数
	int *row_counts = (int *) calloc(nrow, sizeof(int));
	int *col_counts = (int *) calloc(ncolA, sizeof(int));

	int i, j, k;
	int max_row1s = 0;					// 记录初始矩阵里行中非零元素数目的最大值
	int max_col1s = 0;					// 记录初始矩阵里列中非零元素数目的最大值
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncolA; j++) {
			if (A[i][j] != 0) {
				row_counts[i] += 1;
				if (row_counts[i] > max_row1s)
					max_row1s = row_counts[i];

				col_counts[j] += 1;
				if (col_counts[j] > max_col1s)
					max_col1s = col_counts[j];
			}
		}
	}
	//printf("max_row1s: %d, max_col1s: %d\n", max_row1s, max_col1s);
	//for (i=0; i<ncolA; i++)
	//	printf("%d\t", col_counts[i]);
	//printf("\n");

	//ssList **RowID_lists;
	//ssList **ColID_lists;
	//RowID_lists = (ssList *) malloc(sizeof(ssList*) * (max_row1s+1));
	ssList **RowID_lists = (ssList **) malloc(sizeof(ssList*) * (max_row1s+1)); // 指针数组，每个指针指向一个双向链表，同一个链表中的行具有相同数目的非零元素
	ssList **ColID_lists = (ssList **) malloc(sizeof(ssList*) * (max_col1s+1)); // 指针数组，每个指针指向一个双向链表，同一个链表中的列具有相同数目的非零元素
	for (i=0; i<max_row1s+1; i++)
	{
		RowID_lists[i] = (ssList *) malloc(sizeof(ssList));
		RowID_lists[i]->ssFirst = RowID_lists[i]->ssLast = NULL;
	}
	for (i=0; i<max_col1s+1; i++)
	{
		ColID_lists[i] = (ssList *) malloc(sizeof(ssList));
		ColID_lists[i]->ssFirst = ColID_lists[i]->ssLast = NULL;
	}
	// 分类具有不同数目非零元素的行和列
	//
	int allzero_rows = 0;
	for (i=0; i<nrow; i++)
	{
		int row_count = row_counts[i];
		Subscript *rowID = malloc(sizeof(Subscript));
		rowID->index = i;
		rowID->nonzeros = row_count;
		rowID->prev = rowID->next = NULL;
		insertSubAtBeginning(RowID_lists[row_count], rowID);
		if (row_count == 0)
			allzero_rows += 1;
	}
	int allzero_cols = 0;
	for (i=0; i<ncolA; i++)
	{
		int col_count = col_counts[i];
		Subscript *colID = malloc(sizeof(Subscript));
		colID->index = i;
		colID->nonzeros = col_count;
		colID->prev = colID->next = NULL;
		insertSubAtBeginning(ColID_lists[col_count], colID);
		if (col_count == 0)
			allzero_cols += 1;
	}

	int nonallzero_rows = nrow - allzero_rows;
	int nonallzero_cols = ncolA - allzero_cols;
	printf("there are %d all-zero rows and %d all-zero cols in the matrix.\n", allzero_rows, allzero_cols);
	//printf("there are %d not all-zero rows and %d not all-zero cols in the matrix.\n", nonallzero_rows, nonallzero_cols);
	//Jan28注：如果有全零列，则必须在pivoting最后把它们的序号加上去

	// 有了row_counts[], col_counts[], RowID_lists, ColID_lists之后，即开始Markowitz reordering
	int pivots_found = 0;
	int removed_cols = 0;
	int toallzero_cols = 0;
	int removed_rows = 0;
	int toallzero_rows = 0;
	int singletons = 0;
	while (pivots_found != ncolA)
	{
		//printf("%d rows removed, %d reduced to all-zero.\n", removed_rows, toallzero_rows);
		//printf("%d cols removed, %d reduced to all-zero.\n", removed_cols, toallzero_cols);

		// 计算Markowitz count
		int potential_r  = -1;
		int potential_c  = -1;
		int potential_mc = -1;				// potential Markowitz count
		int current_mc;
		int mc_minipos;				// the minimum possible Markowitz count in each row test
		int direct_found = 0;
		// search for one pivot
		Subscript *rows_inchecking, *cols_inchecking;
		int row_id, col_id;
		
		// 按行搜索，并且逐次遍历
		int searched_rows = 0;
		for (i=1; i<=max_row1s; i++)
		{
			// 先按行的非零元素多少搜索
			mc_minipos = (i - 1) * (i - 1);
			if (RowID_lists[i]->ssFirst == NULL)
				continue;

			rows_inchecking = RowID_lists[i]->ssFirst;
			while (rows_inchecking != NULL && (searched_rows < ZLATEVS))
			{
				searched_rows += 1;
				row_id = rows_inchecking->index;
				//printf("row %d contains only one entry.\n", row_id);

				// 再按列的非零元素多少做test
				for (j=1; j<=max_col1s; j++)
				{
					cols_inchecking = ColID_lists[j]->ssFirst;
					while (cols_inchecking != NULL)
					{
						col_id = cols_inchecking->index;
						if (A[row_id][col_id] != 0)
						{
							// we have found an entry in (row_id)-th row
							current_mc = (i-1) * (j-1);
							if (current_mc == 0)
							{
								
								potential_r = row_id;
								potential_c = col_id;
								if (i==1)
									singletons += 1;
								//	printf("a singleton row is found.\n");
								goto found;
							}
							else if (potential_mc == -1 || current_mc < potential_mc)
							{
								potential_r = row_id;
								potential_c = col_id;
								potential_mc = current_mc;
							}
						}
						cols_inchecking = cols_inchecking->next;
					}
				}
				rows_inchecking = rows_inchecking->next;
			}
		}

		
		if (potential_r == -1 || potential_c == -1)
		{
			//printf("error: no pivot is found, the matrix is not full-rank.\n");
			printf("%d rows reduced to all-zero.\n", toallzero_rows);
			printf("%d cols reduced to all-zero.\n", toallzero_cols);
			printf("(partial success) %d/%d pivots were found out of %d rows.\n", pivots_found, ncolA, nrow);
			// 说明有被reduce到全零的行列，把他们当作pivot，anyway
			Subscript *sub_ptt;
			sub_ptt = ColID_lists[0]->ssFirst;
			int zerocols = 0;
			while(sub_ptt != NULL)
			{
				zerocols += 1;
				Subscript *newCpivot0 = malloc(sizeof(Subscript));
				newCpivot0->index = sub_ptt->index;
				newCpivot0->nonzeros = sub_ptt->nonzeros;
				newCpivot0->prev = newCpivot0->next = NULL;
				insertSubAtEnd(ColPivots, newCpivot0);
				sub_ptt = sub_ptt->next;
			}

			sub_ptt = RowID_lists[0]->ssFirst;
			//for (i=0; i<(toallzero_rows+allzero_cols); i++)
			for (i=0; i<zerocols; i++)
			{
				Subscript *newRpivot0 = malloc(sizeof(Subscript));
				newRpivot0->index = sub_ptt->index;
				newRpivot0->nonzeros = sub_ptt->nonzeros;
				newRpivot0->prev = newRpivot0->next = NULL;
				insertSubAtEnd(RowPivots, newRpivot0);
				sub_ptt = sub_ptt->next;
			}

			printf("There are %d/%d singleton rows were found as pivots.\n", singletons, pivots_found);
			//return (pivots_found+toallzero_cols);
			return ncolA;
		}

found:	
		// 找到了一个pivot，保存之
		//printf("a pivot at (%d, %d) is found.\n", potential_r, potential_c);
		pivots_found += 1;
		int p_r = potential_r;
		int p_c = potential_c;
		
		// 保存该pivot的坐标
		Subscript *newRpivot = malloc(sizeof(Subscript));
		newRpivot->index = p_r;
		newRpivot->nonzeros = row_counts[p_r];
		newRpivot->prev = newRpivot->next = NULL;
		insertSubAtEnd(RowPivots, newRpivot);
		Subscript *newCpivot = malloc(sizeof(Subscript));
		newCpivot->index = p_c;
		newCpivot->nonzeros = col_counts[p_c];
		newCpivot->prev = newCpivot->next = NULL;
		insertSubAtEnd(ColPivots, newCpivot);

		// 更新row_counts[], col_counts[], RowID_lists, ColID_lists
		// 1, check potential_r 这一行的所有非零元素，更改它们对应的列的非零元素数目(col_counts[]以及调整ColID_lists)
		int nzs;
		Subscript *ss_pt;
		Subscript *ss_pt_next;
		// 注意有些行和列已经被eliminated了，因此这里采取遍历Subscript对象来更新
		// 1, Update columns
		for (i=1; i<=max_col1s; i++)
		{
			ss_pt = ColID_lists[i]->ssFirst;
			while (ss_pt != NULL)
			{
				if (A[p_r][ss_pt->index] != 0)
				{
					// 该subscript对象将要被更新处理，因此要记下当前遍历的位置
					ss_pt_next = ss_pt->next;
					if (ss_pt->index == p_c)
					{
						removeSubscript(ColID_lists[i], ss_pt);
						removed_cols += 1;
						col_counts[ss_pt->index] -= 1;
						free(ss_pt);
					}
					else
					{
						removeSubscript(ColID_lists[i], ss_pt);
						ss_pt->nonzeros -= 1;
						insertSubAtBeginning(ColID_lists[ss_pt->nonzeros], ss_pt);
						if (ss_pt->nonzeros == 0)
							toallzero_cols += 1;
						col_counts[ss_pt->index] -= 1;
					}
					ss_pt = ss_pt_next;
				}
				else
				{
					ss_pt = ss_pt->next;
				}
			}
		}
		// 2, Update rows
		for (j=1; j<=max_row1s; j++)
		{
			ss_pt = RowID_lists[j]->ssFirst;
			while (ss_pt != NULL)
			{
				if (A[ss_pt->index][p_c] != 0)
				{
					// 该subscript对象将要被更新处理，因此要记下当前遍历的位置
					ss_pt_next = ss_pt->next;
					if (ss_pt->index == p_r)
					{
						removeSubscript(RowID_lists[j], ss_pt);
						removed_rows += 1;
						row_counts[ss_pt->index] -= 1;
						free(ss_pt);
					}
					else
					{
						removeSubscript(RowID_lists[j], ss_pt);
						ss_pt->nonzeros -= 1;
						if (ss_pt->nonzeros == 0)
							toallzero_rows += 1;
						insertSubAtBeginning(RowID_lists[ss_pt->nonzeros], ss_pt);
						row_counts[ss_pt->index] -= 1;
					}
					ss_pt = ss_pt_next;
				}
				else
				{
					ss_pt = ss_pt->next;
				}
			}
		}

	}

	free(row_counts);
	free(col_counts);
	for (i=0; i<max_row1s+1; i++)
		free_subscriptList(RowID_lists[i]);
	for (i=0; i<max_col1s+1; i++)
		free_subscriptList(ColID_lists[i]);

	printf("%d rows reduced to all-zero.\n", toallzero_rows);
	printf("%d cols reduced to all-zero.\n", toallzero_cols);
	printf("(full success) %d/%d pivots were found out of %d rows.\n", pivots_found, ncolA, nrow);
	printf("There are %d/%d singleton rows were found as pivots.\n", singletons, pivots_found);
	return pivots_found;
}

// use progressive inactivation to do pivoting:
//  1) inactivate some columns and only perform pivoting on the rest of the "active" sub-matrix
//  2) if singleton row cannot be found in the middle of pivoting, declare more inactive columns
//  3) given the structure (heavier columns are in the back), declare inactive columns from the back
static int Inactivation_pivoting(int nrow, int ncolA, GF_ELEMENT *A[], ssList *RowPivots, ssList *ColPivots)
{
	printf("Pivoting matrix of size %d x %d via inactivation.\n", nrow, ncolA);

	int i, j, k;
	// 对矩阵中初始非零元素进行计数
	int *row_counts;
	int *col_counts;
	row_counts = (int *) malloc(nrow  * sizeof(int));
	col_counts = (int *) malloc(ncolA * sizeof(int));
	memset(row_counts, 0, nrow  * sizeof(int));
	memset(col_counts, 0, ncolA * sizeof(int));

	int max_row1s = 0;					// 记录初始矩阵里行中非零元素数目的最大值
	int max_col1s = 0;					// 记录初始矩阵里列中非零元素数目的最大值
	for (i=0; i<nrow; i++)
	{
		for (j=0; j<ncolA; j++)
		{
			if (A[i][j] != 0)
			{
				row_counts[i] += 1;
				if (row_counts[i] > max_row1s)
					max_row1s = row_counts[i];

				col_counts[j] += 1;
				if (col_counts[j] > max_col1s)
					max_col1s = col_counts[j];
			}
		}
	}
	// 创建指针数组，每个指针指向一个双向链表，同一个链表中的列具有相同数目的非零元素
	// 该链表用来保存active的列的标号，按非零元素个数排列是为了方便inactivate含非零元素最多的列
	ssList **ColID_lists = (ssList **) malloc(sizeof(ssList*) * (max_col1s+1)); 
	for (i=0; i<max_col1s+1; i++)
	{
		ColID_lists[i] = (ssList *) malloc(sizeof(ssList));
		ColID_lists[i]->ssFirst = ColID_lists[i]->ssLast = NULL;
	}
	for (i=0; i<ncolA; i++)
	{
		int col_count = col_counts[i];
		Subscript *colID = malloc(sizeof(Subscript));
		colID->index = i;
		colID->nonzeros = col_count;
		colID->prev = colID->next = NULL;
		insertSubAtBeginning(ColID_lists[col_count], colID);
	}

	// three possible "states" of each column
	// 0 - active
	// 1 - inactivated
	// 2 - removed because an entry of it was chosen as the pivot
	GF_ELEMENT *inactives;				// an array record the state of columns
	inactives = (GF_ELEMENT *) malloc(sizeof(GF_ELEMENT) * ncolA);
	memset(inactives, 0, sizeof(GF_ELEMENT)*ncolA);

	int inactivated = 0;		// number of inactivated columns
	int active = ncolA;			// number of active columns
	
	int p_r;					// used to store row index of the chosen pivot 
	int p_c;					// used to store col index of the chosen pivot
	int singleton_r_found;
	int ia;
	int selected_pivots = 0;
	while (active != 0)
	{
		p_r = -1;
		p_c = -1;
		singleton_r_found = 0;
		for (i=0; i<nrow; i++)
		{
			if (row_counts[i] == 1)
			{
				singleton_r_found = 1;
				//printf("singleton row is found at row %d.\n", i);
				p_r = i;
				break;
			}
		}

		if (singleton_r_found == 1)
		{
			// a singleton row is found, store the pivot
			for (j=0; j<ncolA; j++)
			{
				if ((inactives[j]==0) && (A[p_r][j]!=0))
				{
					p_c = j;
					break;
				}
			}
			if (p_c == -1)
			{
				printf("error: failed to find the nonzero element in the singlton row.\n");
				exit(1);
			}
			// save the pivot
			// 保存该pivot的坐标
			Subscript *newRpivot = malloc(sizeof(Subscript));
			newRpivot->index = p_r;
			newRpivot->nonzeros = 1;
			newRpivot->prev = newRpivot->next = NULL;
			insertSubAtEnd(RowPivots, newRpivot);
			Subscript *newCpivot = malloc(sizeof(Subscript));
			newCpivot->index = p_c;
			newCpivot->nonzeros = col_counts[p_c];
			newCpivot->prev = newCpivot->next = NULL;
			insertSubAtEnd(ColPivots, newCpivot);

			// 更新row_counts
			row_counts[p_r] = -1;			// use -1 to indicate the row has an elelemnt was chosen as pivot
			for (i=0; i<nrow; i++)
			{
				if ( (row_counts[i] != -1) && (A[i][p_c] != 0) )
				{
					row_counts[i] -= 1;
				}
			}
			// 更新ColID_list
			Subscript *ss_pt = ColID_lists[col_counts[p_c]]->ssFirst;
			while (ss_pt != NULL)
			{
				if (ss_pt->index == p_c)
					break;
				ss_pt = ss_pt->next;
			}
			if (ss_pt->index == p_c)
			{
				removeSubscript(ColID_lists[col_counts[p_c]], ss_pt);
				free(ss_pt);
				ss_pt = NULL;
			}

			inactives[p_c] = 2;
			active -= 1;
			selected_pivots += 1;
			//printf("Inactivated: %d Pivots Found: %d\n", inactivated, selected_pivots);
		}
		else
		{
			// no singleton row can be found, declare one column with the most nonzeros as inactive
			// Note: other algorithms may be used to choose a column to inactivate
			for (i=max_col1s; i>=0; i--)
			{
				Subscript *ss_pt = ColID_lists[i]->ssFirst;
				//if ((ss_pt != NULL) && (inactives[ss_pt->index] == 0))
				if (ss_pt != NULL)
				{
					inactives[ss_pt->index] = 1;
					inactivated += 1;
					active -= 1;
					for (k=0; k<nrow; k++)
					{
						if ((A[k][ss_pt->index] != 0) && (row_counts[k] != -1))
							row_counts[k] -= 1;
					}

					removeSubscript(ColID_lists[i], ss_pt);
					free(ss_pt);
					ss_pt = NULL;
					break;
				}
			}
			selected_pivots += 1;
			//printf("Inactivated: %d Pivots Found: %d\n", inactivated, selected_pivots);
		}

	}
	// assign pivots for the dense part (any ordering is fine)
	for (i=0; i<ncolA; i++)
	{
		if (inactives[i] == 1)
		{
			// find an arbitray row that has not been selected 
			int j_candidate;
			for (j=0; j<nrow; j++)
			{
				if (row_counts[j] != -1)
				{
					j_candidate = j;
					if (A[j][i] != 0)
						break;
				}
			}
			j = j_candidate;

			// 保存该pivot的坐标
			Subscript *newRpivot = malloc(sizeof(Subscript));
			newRpivot->index = j;
			newRpivot->nonzeros = 0;
			newRpivot->prev = newRpivot->next = NULL;
			insertSubAtEnd(RowPivots, newRpivot);
			Subscript *newCpivot = malloc(sizeof(Subscript));
			newCpivot->index = i;
			newCpivot->nonzeros = 0;
			newCpivot->prev = newCpivot->next = NULL;
			insertSubAtEnd(ColPivots, newCpivot);
			row_counts[j] = -1;
			inactives[i] = 2;
			//printf("pivot at (%d, %d) is chosen for the dense part.\n", j, i);
		}
	}

	// free up memories
	free(inactives);

	for (i=0; i<max_col1s+1; i++)
		free_subscriptList(ColID_lists[i]);
	free(ColID_lists);

	free(row_counts);
	free(col_counts);
	return inactivated;
}

// Given pivot sequence, re-order matrices
static void matrices_reordering(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[], ssList *RowPivots, ssList *ColPivots, int npv)
{
	int i, j, k;
	int row_id, col_id;
	Subscript *sub_row, *sub_col;
	Subscript *ss_pt;
	sub_row = RowPivots->ssFirst;
	sub_col = ColPivots->ssFirst;
	GF_ELEMENT temp;
	for (k=0; k<npv; k++)
	{
		row_id = sub_row->index;
		col_id = sub_col->index;

		// swap rows first
		if (row_id != k)
		{
			for (j=0; j<ncolA; j++)
			{
				temp = A[k][j];
				A[k][j] = A[row_id][j];
				A[row_id][j] = temp;
			}
			// For row operations, B needs to be swapped as well
			for (j=0; j<ncolB; j++)
			{
				temp = B[k][j];
				B[k][j] = B[row_id][j];
				B[row_id][j] = temp;
			}

			// due to the swap of rows, update RowPivots to track positions of the rest pivots
			ss_pt = sub_row->next;
			while (ss_pt != NULL)
			{
				if (ss_pt->index == k)
					ss_pt->index = row_id;

				ss_pt = ss_pt->next;
			}
		}

		// swap columns next
		if (col_id != k)
		{
			for (i=0; i<nrow; i++)
			{
				temp = A[i][k];
				A[i][k] = A[i][col_id];
				A[i][col_id] = temp;
			}
			// due to the swap of columns, update ColPivots to track positions of the rest pivots
			ss_pt = sub_col->next;
			while (ss_pt != NULL)
			{
				if (ss_pt->index == k)
					ss_pt->index = col_id;

				ss_pt = ss_pt->next;
			}
		}

		sub_row = sub_row->next;
		sub_col = sub_col->next;
	}
}

// Re-order columns of the matrix
static void permute_matrix_columns(int nrow, int ncolA, GF_ELEMENT *A[], ssList *ColPivots)
{
	int i, j, k;
	int col_id;
	Subscript *sub_col;
	Subscript *ss_pt;
	sub_col = ColPivots->ssFirst;
	GF_ELEMENT temp;
	for (k=0; k<ncolA; k++)
	{
		col_id = sub_col->index;

		// swap columns
		if (col_id != k)
		{
			for (i=0; i<nrow; i++)
			{
				temp = A[i][k];
				A[i][k] = A[i][col_id];
				A[i][col_id] = temp;
			}
			// due to the swap of columns, update ColPivots to track positions of the rest pivots
			ss_pt = sub_col->next;
			while (ss_pt != NULL)
			{
				if (ss_pt->index == k)
					ss_pt->index = col_id;

				ss_pt = ss_pt->next;
			}
		}
		sub_col = sub_col->next;
	}
}

// perform back-substitution on full-rank upper trianguler matrix A
/*
static long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[])
{
	//printf("entering back_substitute()...\n");
	long long operations = 0;

	// 1, transform the upper triangular matrix A into diagonal.
	int i, j, k, l;
	for (i=ncolA-1; i>=0; i--)
	{
		// eliminate all items above A[i][i]
		#pragma omp parallel for
		for (j=0; j<i; j++)
		{
			if (A[j][i] == 0)
				continue;				// the matrix we are dealing here is high likely being sparse, so we use this to avoid unnecessary operations
			GF_ELEMENT quotient = galois_divide(A[j][i], A[i][i], GF_ORDER);
			operations += 1;

			A[j][i] = galois_sub(A[j][i], galois_multiply(A[i][i], quotient, GF_ORDER));
			operations += 1;

			// doing accordingly to B
			galois_multiply_add_region(B[i], B[j], quotient, ncolB, GF_ORDER);
			operations += ncolB;
			//operations += 1;
		}
	}

	// 2, transform matrix A into identity matrix
	for (l=0; l<ncolA; l++)
	{
		if (A[l][l] != 0 && A[l][l] != 1)
		{
			for (k=0; k<ncolB; k++)
			{
				B[l][k] = galois_divide(B[l][k], A[l][l], GF_ORDER);
				operations += 1;
			}
			//operations += 1;
			A[l][l] = 1;
		}
	}

	return operations;
}

// perform forward substitution on a matrix to transform it to a upper triangular structure
static long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[])
{
	//printf("entering foward_substitute()...\n");
	long long operations = 0;
	int i, j, k, m, n, p;
	int pivot;
	GF_ELEMENT quotient;

	// transform A into upper triangular structure by row operation
	int boundary = nrow >= ncolA ? ncolA : nrow;

	for (i=0; i<boundary; i++)
	{
		int has_a_dimension = 1;			// to indicate if this column is all-zeros

		// if the pivot is zero, find a pivot from elements below the diagonal element and swap rows
		if (A[i][i] == 0)
		{
			has_a_dimension = 0;
			if (i == nrow-1)
				break;
			else
			{
				for (pivot=i+1; pivot<nrow; pivot++)
				{
					if (A[pivot][i] != 0)
					{
						has_a_dimension = 1;
						break;
					}
				}
			}
			// if this column is an all-zeros column, just skip this column
			if (has_a_dimension == 0)
			{
				continue;
			}
			else
			{
				// do the swap
				GF_ELEMENT tmp2;
				for (m=0; m<ncolA; m++)
				{
					tmp2 = A[i][m];
					A[i][m] = A[pivot][m];
					A[pivot][m] = tmp2;
				}

				// swap B accordingly
				for (m=0; m<ncolB; m++)
				{
					tmp2 = B[i][m];
					B[i][m] = B[pivot][m];
					B[pivot][m] = tmp2;
				}
			}
		}
		// perform elimination
		for (j=i+1; j<nrow; j++)
		{
			if (A[j][i] == 0)
				continue;			// the matrix we are dealing here is high likely being sparse, so we use this to avoid unnecessary operations
			quotient = galois_divide(A[j][i], A[i][i], GF_ORDER);
			operations += 1;
			// eliminate the items under row i at col i
			galois_multiply_add_region(&(A[i][i]), &(A[j][i]), quotient, ncolA-i, GF_ORDER);
			operations += (ncolA-i);						// Jan27注:此时已经进行到第i列							
			// simultaneously do the same thing on right matrix B
			galois_multiply_add_region(B[i], B[j], quotient, ncolB, GF_ORDER);
			operations += ncolB;
		}
	}
	return operations;
}
*/


// insert a subscript object at the beginning of the list
static void insertSubAtBeginning(ssList *sub_list, Subscript *sub)
{
	if (sub_list->ssFirst == NULL)
	{
		sub_list->ssFirst = sub_list->ssLast = sub;
		sub->prev = NULL;
		sub->next = NULL;
	}
	else
	{
		sub_list->ssFirst->prev = sub;
		sub->next = sub_list->ssFirst;
		sub_list->ssFirst = sub;
		sub->prev = NULL;
	}
}

// insert a subscript object at the end of the list
static void insertSubAtEnd(ssList *sub_list, Subscript *sub)
{
	if (sub_list->ssFirst == NULL && sub_list->ssLast == NULL)
	{
		sub_list->ssFirst = sub_list->ssLast = sub;
		sub->prev = NULL;
		sub->next = NULL;
	}
	else
	{
		sub_list->ssLast->next = sub;
		sub->prev = sub_list->ssLast;
		sub->next = NULL;
		sub_list->ssLast = sub;
	}
}

// remove a subscript object from the list (no free() operation yet)
static void removeSubscript(ssList *sub_list, Subscript *sub)
{
	if (sub->prev == NULL && sub->next == NULL)
	{
		// only one left in the list
		sub_list->ssFirst = sub_list->ssLast = NULL;
	}
	else if (sub->prev == NULL && sub->next != NULL)
	{
		// the element is the head
		sub_list->ssFirst = sub->next;
		sub_list->ssFirst->prev = NULL;
	}
	else if (sub->prev != NULL && sub->next == NULL)
	{
		// the element is the tail
		sub_list->ssLast = sub->prev;
		sub->prev->next = NULL;
	}
	else
	{
		// the element is in the middle of the list
		sub->next->prev = sub->prev;
		sub->prev->next = sub->next;
	}
}

// free the subscript list
static void free_subscriptList(ssList *sub_list)
{
	Subscript *sub = sub_list->ssFirst;
	while (sub != NULL)
	{
		sub_list->ssFirst = sub->next;
		free(sub);
		sub = sub_list->ssFirst;
	}
	free(sub_list);
	sub_list = NULL;
}
