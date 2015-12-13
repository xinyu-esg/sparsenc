/*-----------------------decoderBD.c----------------------
 * Implementation of regular band decoder. It employs pivoting
 * to jointly decode band GNC code and its precode.
 *------------------------------------------------------------*/
#include "common.h"
#include "galois.h"
#include "decoderBD.h"
static int partially_diag_decoding_matrix(struct decoding_context_BD *dec_ctx);
static int apply_parity_check_matrix(struct decoding_context_BD *dec_ctx);
static void finish_recovering_BD(struct decoding_context_BD *dec_ctx);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int *otoc, int *inactives);

// create decoding context for band decoder
void create_dec_context_BD(struct decoding_context_BD *dec_ctx, struct snc_parameter sp)
{
    static char fname[] = "snc_create_dec_context_BD";
    int i, j, k;

    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    if (sp.type != BAND_SNC) {
        fprintf(stderr, "Band decoder only applies to band GNC code.\n");
        return;
    }
    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) 
        fprintf(stderr, "%s: create decoding context failed", fname);

    dec_ctx->sc = sc;

    dec_ctx->finished     = 0;
    dec_ctx->DoF 		  = 0;
    dec_ctx->de_precode   = 0;
    dec_ctx->inactivated  = 0;

    int gensize = dec_ctx->sc->meta.size_g;
    int pktsize = dec_ctx->sc->meta.size_p;
    int numpp   = dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum;

    dec_ctx->coefficient = calloc(numpp, sizeof(GF_ELEMENT*));
    dec_ctx->message     = calloc(numpp, sizeof(GF_ELEMENT*));
    for (i=0; i<numpp; i++) {
        dec_ctx->coefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
        dec_ctx->message[i]     = calloc(pktsize, sizeof(GF_ELEMENT));
    }

    dec_ctx->otoc_mapping = malloc(sizeof(int) * numpp);
    dec_ctx->ctoo_mapping = malloc(sizeof(int) * numpp);
    for (j=0; j<numpp; j++) {
        dec_ctx->otoc_mapping[j]   = j;				// original to current mapping
        dec_ctx->ctoo_mapping[j]   = j;				// current to original mapping
    }

    dec_ctx->overhead     = 0;
    dec_ctx->overheads = calloc(dec_ctx->sc->meta.gnum, sizeof(int));
    dec_ctx->operations   = 0;
}


void process_packet_BD(struct decoding_context_BD *dec_ctx, struct snc_packet *pkt)
{

    dec_ctx->overhead += 1;
    dec_ctx->overheads[pkt->gid] += 1;
    int i, j, k;
    int pivot;
    int pivotfound = 0;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->meta.size_g;
    int pktsize = dec_ctx->sc->meta.size_p;
    int numpp   = dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum;

    // translate GNC encoding vector to full length
    GF_ELEMENT *ces = calloc(numpp, sizeof(GF_ELEMENT));

    if (dec_ctx->de_precode == 0) {
        /*
         * Before precode's check matrix was applied
         */
        for (i=0; i<gensize; i++) {
            int index = dec_ctx->sc->gene[pkt->gid]->pktid[i];
            ces[index] = pkt->coes[i];
        }
        for (i=0; i<numpp; i++) {
            if (ces[i] != 0) {
                if (dec_ctx->coefficient[i][i] != 0) {
                    quotient = galois_divide(ces[i], dec_ctx->coefficient[i][i], GF_POWER);
                    dec_ctx->operations += 1;
                    int band_width = numpp-i > gensize ? gensize : numpp-i;
                    galois_multiply_add_region(ces+i, &(dec_ctx->coefficient[i][i]), quotient, band_width, GF_POWER);
                    dec_ctx->operations += band_width;
                    galois_multiply_add_region(pkt->syms, dec_ctx->message[i], quotient, pktsize, GF_POWER);
                    dec_ctx->operations += pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
    } else {
        /*
         *Parity-check matrix has been applied, and therefore the decoding matrix has been pivoted and re-ordered
         */
        for (i=0; i<gensize; i++) {
            int orig_index = dec_ctx->sc->gene[pkt->gid]->pktid[i];
            int curr_index = dec_ctx->otoc_mapping[orig_index];
            ces[curr_index] = pkt->coes[i];
        }
        for (i=0; i<numpp; i++) {
            if (ces[i] != 0) {
                if (dec_ctx->coefficient[i][i] != 0) {
                    quotient = galois_divide(ces[i], dec_ctx->coefficient[i][i], GF_POWER);
                    dec_ctx->operations += 1;
                    galois_multiply_add_region(ces+i, &(dec_ctx->coefficient[i][i]), quotient, numpp-i, GF_POWER);
                    dec_ctx->operations += (numpp - i);
                    galois_multiply_add_region(pkt->syms, dec_ctx->message[i], quotient, pktsize, GF_POWER);
                    dec_ctx->operations += pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
    }

    //printf("Operations consumed: %lld\n", dec_ctx->operations);

    if (pivotfound == 1) {
        memcpy(dec_ctx->coefficient[pivot], ces, numpp*sizeof(GF_ELEMENT));
        memcpy(dec_ctx->message[pivot], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
        dec_ctx->DoF += 1;
    }
    // If the number of received DoF is equal to NUM_SRC, apply the parity-check matrix.
    // The messages corresponding to rows of parity-check matrix are all-zero.
    if (dec_ctx->DoF == dec_ctx->sc->meta.snum) {
#if defined(GNCTRACE)
        printf("Start to apply the parity-check matrix...\n");
#endif
        int allzeros = partially_diag_decoding_matrix(dec_ctx);	
#if defined(GNCTRACE)
        printf("%d all-zero rows when partially diagonalizing the decoding matrix.\n", allzeros);
#endif
        int missing_DoF = apply_parity_check_matrix(dec_ctx);
#if defined(GNCTRACE)
        printf("After applying the parity-check matrix, %d DoF are missing.\n", missing_DoF);
#endif
        dec_ctx->DoF = numpp - missing_DoF;
        dec_ctx->de_precode = 1;
    }

    if (dec_ctx->DoF == dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum) {
        finish_recovering_BD(dec_ctx);
    }

    snc_free_packet(pkt);
    pkt = NULL;
    free(ces);
    ces = NULL;
}

/* 
 * Partially diagonalize the upper-trianguler decoding matrix, 
 * i.e., remove nonzero elements above nonzero diagonal elements
 */
static int partially_diag_decoding_matrix(struct decoding_context_BD *dec_ctx)
{
    int 		i, j, l;
    int 		nonzero_rows = 0;
    GF_ELEMENT 	quotient;
    long long 	operations = 0;

    int gensize = dec_ctx->sc->meta.size_g;
    int pktsize = dec_ctx->sc->meta.size_p;
    int numpp = dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum;

    int	zero_size	= dec_ctx->sc->meta.cnum + 5;  
    int *zeropivots = (int *) malloc( sizeof(int) * zero_size );	// store indices of columns where the diagonal element is zero
    int zero_p		= 0;										// indicate how many zero pivots have been identified

    for (j=numpp-1; j>=0; j--) {
        if (dec_ctx->coefficient[j][j] == 0) {
            if (zero_p == zero_size) {
                zeropivots = (int *) realloc(zeropivots, sizeof(int) * (zero_size + 5));
                zero_size += 5;
            }
            zeropivots[zero_p++] = j;
            continue;
        } else {
            nonzero_rows += 1;
            int start_row = j-gensize > 0 ? j-gensize : 0;		// the upper triangular form is also in banded form, so no need to go through all rows
            for (i=start_row; i<j; i++) {
                if (dec_ctx->coefficient[i][j] == 0)
                    continue;

                quotient = galois_divide(dec_ctx->coefficient[i][j], dec_ctx->coefficient[j][j], GF_POWER);
                operations += 1;
                dec_ctx->coefficient[i][j] = 0;			// eliminiate the element
                // Important: corresponding operations on behind columns whose diagonal elements are zeros
                for (int z=0; z<zero_p; z++) {
                    l = zeropivots[z];
                    if (dec_ctx->coefficient[j][l] != 0) {
                        dec_ctx->coefficient[i][l] = galois_add(dec_ctx->coefficient[i][l], galois_multiply(dec_ctx->coefficient[j][l], quotient, GF_POWER));
                        operations += 1;
                    }
                }
                // correspoding operations on the message matrix
                galois_multiply_add_region(dec_ctx->message[i], dec_ctx->message[j], quotient, pktsize, GF_POWER);
                operations += pktsize;
            }
        }
    }
    free(zeropivots);
    dec_ctx->operations += operations;
    return (numpp-nonzero_rows);
}

// Apply the parity-check matrix to the decoding matrix; pivot, re-order and try to jointly decode
static int apply_parity_check_matrix(struct decoding_context_BD *dec_ctx)
{
    int i, j, k;
    int num_of_new_DoF = 0;

    int gensize = dec_ctx->sc->meta.size_g;
    int pktsize = dec_ctx->sc->meta.size_p;
    int numpp = dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum;

    // 1, Copy parity-check vectors to the nonzero rows of the decoding matrix
    int p = 0;						// index pointer to the parity-check vector that is to be copyed
    for (i=0; i<numpp; i++) {
        if (dec_ctx->coefficient[i][i] == 0) {
            /* Set the coding vector according to parity-check bits */
            NBR_node *varnode = dec_ctx->sc->graph->l_nbrs_of_r[p]->first;
            while (varnode != NULL) {
                dec_ctx->coefficient[i][varnode->data] = 1;
                varnode = varnode->next;
            }
            dec_ctx->coefficient[i][dec_ctx->sc->meta.snum+p] = 1;
            p++;
            memset(dec_ctx->message[i], 0, sizeof(GF_ELEMENT)*pktsize);			// parity-check vector corresponds to all-zero message
        }
    }

    // 2, Pivot and re-order matrices
    dec_ctx->operations += pivot_matrix_oneround(numpp, numpp, pktsize, dec_ctx->coefficient, dec_ctx->message, dec_ctx->otoc_mapping, &(dec_ctx->inactivated));

    /* Count available innovative rows */
    int missing_DoF = 0;
    for (i=0; i<numpp; i++) {
        if (dec_ctx->coefficient[i][i] == 0)
            missing_DoF++;
        dec_ctx->ctoo_mapping[dec_ctx->otoc_mapping[i]] = i;
    }
    return missing_DoF;
}


// recover decoded packets after NUM_SRC DoF has been received
static void finish_recovering_BD(struct decoding_context_BD *dec_ctx)
{
    int gensize = dec_ctx->sc->meta.size_g;
    int pktsize = dec_ctx->sc->meta.size_p;
    int numpp = dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum;
    //printf("Operations consumed before back substitution: %lld\n", dec_ctx->operations);
    long long bs_ops = back_substitute(numpp, numpp, pktsize, dec_ctx->coefficient, dec_ctx->message);
    dec_ctx->operations += bs_ops;
    //printf("Operations consumed afer back substitution: %lld\n", dec_ctx->operations);
    int i, j, k;
    for (i=0; i<numpp; i++){
        dec_ctx->sc->pp[dec_ctx->ctoo_mapping[i]] = calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dec_ctx->sc->pp[dec_ctx->ctoo_mapping[i]], dec_ctx->message[i], pktsize*sizeof(GF_ELEMENT));
    }
    dec_ctx->finished = 1;
}

void free_dec_context_BD(struct decoding_context_BD *dec_ctx)
{
    for (int i=dec_ctx->sc->meta.snum+dec_ctx->sc->meta.cnum-1; i>=0; i--) {
        free(dec_ctx->coefficient[i]);
        free(dec_ctx->message[i]);
    }
    free(dec_ctx->coefficient);
    free(dec_ctx->message);
    free(dec_ctx->overheads);
    free(dec_ctx->otoc_mapping);
    free(dec_ctx->ctoo_mapping);
    snc_free_enc_context(dec_ctx->sc);
    free(dec_ctx);
    dec_ctx = NULL;
}

