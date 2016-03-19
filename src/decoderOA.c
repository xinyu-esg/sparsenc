/*--------------------------sncOADecoder.c----------------------
 * Implementation of overlap-aware (OA) decoder.
 *-------------------------------------------------------------*/
#include "common.h"
#include "galois.h"
#include "decoderOA.h"

// to store matrices in processing (needed by the decoder)
struct running_matrix
{
    GF_ELEMENT **coefficient;  //[CLASS_SIZE][CLASS_SIZE];
    GF_ELEMENT **message;      //[CLASS_SIZE][EXT_N];
};

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

/* Free running matrix */
static void free_running_matrix(struct running_matrix *mat, int rows);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int *otoc, int *inactives);
extern long pivot_matrix_tworound(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int *otoc, int *inactives);

/*
 * snc_create_dec_context_OA
 * Create context for overlap-aware (OA) decoding
 *  aoh - allowed overhead >=0
 */
struct decoding_context_OA *create_dec_context_OA(struct snc_parameters *sp, int aoh)
{
    static char fname[] = "snc_create_dec_context_OA";
    int i, j, k;

    struct decoding_context_OA *dec_ctx;
    if ((dec_ctx = malloc(sizeof(struct decoding_context_OA))) == NULL)
        return NULL;

    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) {
        fprintf(stderr, "%s: create decoding context failed", fname);
        goto AllocError;
    }

    dec_ctx->sc = sc;

    dec_ctx->aoh        = aoh;
    dec_ctx->finished   = 0;
    dec_ctx->OA_ready   = 0;
    dec_ctx->local_DoF  = 0;
    dec_ctx->global_DoF = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;

    // Allocate matrices for per-generation decoding
    dec_ctx->Matrices = calloc(dec_ctx->sc->params.gnum, sizeof(struct running_matrix*));
    if (dec_ctx->Matrices == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->Matrices\n", fname);
        goto AllocError;
    }
    for (i=0; i<dec_ctx->sc->params.gnum; i++) {
        dec_ctx->Matrices[i] = calloc(1, sizeof(struct running_matrix));
        if (dec_ctx->Matrices[i] == NULL) {
            fprintf(stderr, "%s: malloc dec_ctx->Matrices[%d]\n", fname, i);
            goto AllocError;
        }
        // Allocate coefficient and message matrices in running_matrix
        // coefficeint: size_g x size_g
        // message:     size_g x size_p
        // Dim-1) Pointers to each row
        dec_ctx->Matrices[i]->coefficient = calloc(gensize, sizeof(GF_ELEMENT*));
        if (dec_ctx->Matrices[i]->coefficient == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->coefficient\n", fname, i);
            goto AllocError;
        }
        dec_ctx->Matrices[i]->message = calloc(gensize, sizeof(GF_ELEMENT*));
        if (dec_ctx->Matrices[i]->message == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage\n", fname, i);
            goto AllocError;
        }
        // Dim-2) Elements of each row
        for (j=0; j<gensize; j++) {
            dec_ctx->Matrices[i]->coefficient[j] = calloc(gensize, sizeof(GF_ELEMENT));
            if (dec_ctx->Matrices[i]->coefficient[j] == NULL) {
                fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->coefficient[%d]\n", fname, i, j);
                goto AllocError;
            }
            dec_ctx->Matrices[i]->message[j] = calloc(pktsize, sizeof(GF_ELEMENT));
            if (dec_ctx->Matrices[i]->message[j] == NULL) {
                fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage[%d]\n", fname, i, j);
                goto AllocError;
            }
        }
    }

    /*
     * We don't allocate memory for global decoding (ie GDM) here. We only allocate
     * when OA ready. This avoids occupying a big amount of memory for a long time.
     */

    // performance indices
    dec_ctx->operations = 0;
    dec_ctx->overhead   = 0;
    dec_ctx->ops1 = dec_ctx->ops2 = dec_ctx->ops3 = dec_ctx->ops4 = 0;
    return dec_ctx;

AllocError:
    free_dec_context_OA(dec_ctx);
    dec_ctx = NULL;
    return NULL;
}

void process_packet_OA(struct decoding_context_OA *dec_ctx, struct snc_packet *pkt)
{
    dec_ctx->overhead += 1;

    int i, j, k;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;

    // start processing
    int gid = pkt->gid;
    int pivotfound = 0;
    int pivot;

    /*
     * If decoder is not OA ready, process the packet within the generation.
     */
    if (dec_ctx->OA_ready != 1) {
        GF_ELEMENT *pkt_coes = calloc(gensize, sizeof(GF_ELEMENT));
        if (dec_ctx->sc->params.bnc) {
            for (i=0; i<gensize; i++)
                pkt_coes[i] = get_bit_in_array(pkt->coes, i);
        } else {
            memcpy(pkt_coes, pkt->coes, gensize*sizeof(GF_ELEMENT));
        }
        // Translate the encoding vector to the sorted form as in the generation
        struct running_matrix *matrix = dec_ctx->Matrices[gid];
        for (i=0; i<gensize; i++) {
            if (pkt_coes[i] != 0) {
                if (matrix->coefficient[i][i] != 0) {
                    quotient = galois_divide(pkt_coes[i], matrix->coefficient[i][i], GF_POWER);
                    galois_multiply_add_region(&(pkt_coes[i]), &(matrix->coefficient[i][i]), quotient, gensize-i, GF_POWER);
                    galois_multiply_add_region(pkt->syms, matrix->message[i], quotient, pktsize, GF_POWER);
                    dec_ctx->operations += 1 + (gensize - i) + pktsize;
                    dec_ctx->ops1 += 1 + (gensize -i) + pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
        // cache as normal GNC packet
        if (pivotfound == 1) {
            memcpy(matrix->coefficient[pivot], pkt_coes, gensize*sizeof(GF_ELEMENT));
            memcpy(matrix->message[pivot], pkt->syms, pktsize*sizeof(GF_ELEMENT));
            dec_ctx->local_DoF += 1;
        }
        free(pkt_coes);

        if ((dec_ctx->local_DoF >= dec_ctx->sc->params.snum) && (dec_ctx->overhead >= (dec_ctx->sc->params.snum+dec_ctx->aoh))) {
            dec_ctx->OA_ready = 1;
            // When OA ready, convert LDMs to upper triangular form
            long ops = running_matrix_to_REF(dec_ctx);
            dec_ctx->operations += ops;
            dec_ctx->ops1 += ops;
            // Combine LDMs to GDM and apply inactivation pivoting
            construct_GDM(dec_ctx);

            // If numpp innovative packets are received, recover all
            // source packets from JMBcoeffcient and JMBmessage
            if (dec_ctx->global_DoF == numpp)
                diagonalize_GDM(dec_ctx);
        }
    } else {
        /*
         * If decoder is OA ready, process the packet against the global matrix directly.
         * Since the decoder is OA ready, translate the local encoded vector (LEV)
         * to global encoding vector (GEV). Since the GDM was probably pivoted, need
         * to transform the GEV according to the pivoting order.
         */
        GF_ELEMENT *re_ordered = calloc(numpp, sizeof(GF_ELEMENT));
        for (i=0; i<gensize; i++) {
            /* obtain current index position of pktid */
            int curr_pos = dec_ctx->otoc_mapping[dec_ctx->sc->gene[gid]->pktid[i]];
            if (dec_ctx->sc->params.bnc) {
                re_ordered[curr_pos] = get_bit_in_array(pkt->coes, i);
            } else {
                re_ordered[curr_pos] = pkt->coes[i];
            }
        }

        /*
         * Process the reordered GEV against GDM
         */
        pivot = -1;
        for (int m=0; m<numpp; m++) {
            if (re_ordered[m] != 0) {
                if (dec_ctx->JMBcoefficient[m][m] != 0) {
                    // mask the encoding vector and message over the JMB decoding matrix
                    GF_ELEMENT quotient = galois_divide(re_ordered[m], dec_ctx->JMBcoefficient[m][m], GF_POWER);
                    dec_ctx->operations += 1;
                    dec_ctx->ops3 += 1;
                    if ( m < (numpp - dec_ctx->inactives) ) {
                        // Only needs to be multiply-and-add to the inactive part
                        // this saves computation
                        int maa_start = numpp - dec_ctx->inactives;
                        galois_multiply_add_region(re_ordered+maa_start, &(dec_ctx->JMBcoefficient[m][maa_start]), quotient, dec_ctx->inactives, GF_POWER);
                        dec_ctx->operations += dec_ctx->inactives;
                        dec_ctx->ops3 += dec_ctx->inactives;
                        re_ordered[m] = 0;
                    } else {
                        galois_multiply_add_region(re_ordered+m, &(dec_ctx->JMBcoefficient[m][m]), quotient, numpp-m, GF_POWER);
                        dec_ctx->operations += (numpp - m);
                        dec_ctx->ops3 += (numpp - m);
                    }
                    galois_multiply_add_region(pkt->syms, dec_ctx->JMBmessage[m], quotient, pktsize, GF_POWER);
                    dec_ctx->operations += pktsize;
                    dec_ctx->ops3 += pktsize;
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
    if (dec_ctx->finished && get_loglevel() == TRACE) {
        printf("OA splitted operations: %.2f %.2f %.2f %.2f\n",
                (double) dec_ctx->ops1/dec_ctx->sc->params.snum/pktsize,
                (double) dec_ctx->ops2/dec_ctx->sc->params.snum/pktsize,
                (double) dec_ctx->ops3/dec_ctx->sc->params.snum/pktsize,
                (double) dec_ctx->ops4/dec_ctx->sc->params.snum/pktsize);
    }
}

void free_dec_context_OA(struct decoding_context_OA *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
    if (dec_ctx->sc != NULL)
        snc_free_enc_context(dec_ctx->sc);
    int i, j, k;
    if (dec_ctx->Matrices != NULL) {
        for (i=0; i<dec_ctx->sc->params.gnum; i++){
            // Free each decoding matrix
            if (dec_ctx->Matrices[i] != NULL)
                free_running_matrix(dec_ctx->Matrices[i], dec_ctx->sc->params.size_g);
            dec_ctx->Matrices[i] = NULL;
        }
        free(dec_ctx->Matrices);
        dec_ctx->Matrices = NULL;
    }
    if (dec_ctx->JMBcoefficient != NULL) {
        for (j=0; j<dec_ctx->sc->params.snum+dec_ctx->sc->params.cnum+dec_ctx->aoh; j++) {
            if (dec_ctx->JMBcoefficient[j] != NULL)
                free(dec_ctx->JMBcoefficient[j]);
        }
        free(dec_ctx->JMBcoefficient);
    }
    if (dec_ctx->JMBmessage != NULL) {
        for (j=0; j<dec_ctx->sc->params.snum+dec_ctx->sc->params.cnum+dec_ctx->aoh; j++) {
            if (dec_ctx->JMBmessage[j] != NULL)
                free(dec_ctx->JMBmessage[j]);
        }
        free(dec_ctx->JMBmessage);
    }
    if (dec_ctx->otoc_mapping != NULL)
        free(dec_ctx->otoc_mapping);
    if (dec_ctx->ctoo_mapping != NULL)
        free(dec_ctx->ctoo_mapping);
    free(dec_ctx);
    dec_ctx = NULL;
    return;
}

static void free_running_matrix(struct running_matrix *mat, int rows)
{
    if (mat != NULL) {
        if (mat->coefficient != NULL) {
            for (int i=0; i<rows; i++) {
                if (mat->coefficient[i] != NULL)
                    free(mat->coefficient[i]);
            }
            free(mat->coefficient);
        }
        mat->coefficient = NULL;
        if (mat->message != NULL) {
            for (int i=0; i<rows; i++) {
                if (mat->message[i] != NULL)
                    free(mat->message[i]);
            }
            free(mat->message);
        }
        mat->message = NULL;
        return;
    }
}

static void diagonalize_GDM(struct decoding_context_OA *dec_ctx)
{
    static char fname[] = "finish_recovering_inactivation";
    int i, j, k;
    int pos;
    int pktid;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;

    // Recover inactivated packets
    if (get_loglevel() == TRACE) {
        printf("Finishing decoding...\n");
        printf("Recovering \"inactive\" packets...\n");
    }
    int ias = dec_ctx->inactives;
    GF_ELEMENT **ces_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ias; i++) {
        ces_submatrix[i] = calloc(ias, sizeof(GF_ELEMENT));
        memcpy(ces_submatrix[i], &(dec_ctx->JMBcoefficient[numpp-ias+i][numpp-ias]), ias*sizeof(GF_ELEMENT));
    }

    /* Perform back substitution to reduce the "ias x ias" matrix to identity matrix */
    long long ops = back_substitute(ias, ias, pktsize, ces_submatrix, &(dec_ctx->JMBmessage[numpp-ias]));
    dec_ctx->operations += ops;
    dec_ctx->ops4 += ops;

    // Recover decoded overlapping packets
    for (i=0; i<ias; i++) {
        // get original pktid at column (numpp-ias+i0
        pktid = dec_ctx->ctoo_mapping[numpp-ias+i];
        // Construct decoded packets
        if ( (dec_ctx->sc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
            fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, pktid);
        memcpy(dec_ctx->sc->pp[pktid], dec_ctx->JMBmessage[numpp-ias+i], sizeof(GF_ELEMENT)*pktsize);
    }
    /* free ces_submatrix */
    for (i=0; i<ias; i++)
        free(ces_submatrix[i]);
    free(ces_submatrix);


    // Recover active packets
    if (get_loglevel() == TRACE)
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
                galois_multiply_add_region(dec_ctx->JMBmessage[i], dec_ctx->sc->pp[pktid], quotient, pktsize, GF_POWER);
                dec_ctx->JMBcoefficient[i][j] = 0;
                dec_ctx->operations += pktsize;
                dec_ctx->ops4 += pktsize;
            }
        }

        // Convert diagonal elements of top-left part of T to 1
        quotient = dec_ctx->JMBcoefficient[i][i];
        if (quotient != 1) {
            galois_multiply_region(dec_ctx->JMBmessage[i], galois_divide(1, quotient, GF_POWER), pktsize, GF_POWER);
            dec_ctx->operations += pktsize;
            dec_ctx->ops4 += pktsize;
            dec_ctx->JMBcoefficient[i][i] = 1;
        }

        // Save the decoded packet
        pktid = dec_ctx->ctoo_mapping[i];
        if ( dec_ctx->sc->pp[pktid] != NULL )
            fprintf(stderr, "%s：warning: packet %d is already recovered.\n", fname, pktid);
        if ( (dec_ctx->sc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
            fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, pktid);
        memcpy(dec_ctx->sc->pp[pktid], dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT)*pktsize);
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

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;


    for (i=0; i<dec_ctx->sc->params.gnum; i++) {
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

                quotient = galois_divide(matrix->coefficient[l][k], matrix->coefficient[k][k], GF_POWER);
                operations += 1;
                matrix->coefficient[l][k] = 0;
                // 注意后面有的列可能并非为零列(也就是对角元素非零)，这些也要eliminate
                for (int m=k+1; m<gensize; m++) {
                    if (matrix->coefficient[m][m] == 0) {
                        matrix->coefficient[l][m] = galois_add(matrix->coefficient[l][m], galois_multiply(matrix->coefficient[k][m], quotient, GF_POWER));
                        operations += 1;
                    }
                }
                galois_multiply_add_region(matrix->message[l], matrix->message[k], quotient, pktsize, GF_POWER);
                operations += pktsize;
            }
        }
        if (get_loglevel() == TRACE) {
            if (consecutive == 1)
                printf("Class %d is self-decodable.\n", i);
        }
    }
    return operations;
}

static void construct_GDM(struct decoding_context_OA *dec_ctx)
{
    static char fname[] = "construct_GDM";
    int i, j, k;
    struct running_matrix *matrix;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;

    //Allocate GDM to snc_dec_context, apply precoding matrix
    dec_ctx->JMBcoefficient = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
    if (dec_ctx->JMBcoefficient == NULL)
        fprintf(stderr, "%s: calloc dec_ctx->JMBcoefficient\n", fname);
    dec_ctx->JMBmessage     = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
    if (dec_ctx->JMBmessage == NULL)
        fprintf(stderr, "%s: calloc dec_ctx->JMBmessage\n", fname);

    for (i=0; i<numpp+dec_ctx->aoh; i++) {
        dec_ctx->JMBcoefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
        dec_ctx->JMBmessage[i]     = calloc(pktsize, sizeof(GF_ELEMENT));
    }
    dec_ctx->inactives   = 0;
    dec_ctx->otoc_mapping = malloc(sizeof(int) * numpp);
    dec_ctx->ctoo_mapping = malloc(sizeof(int) * numpp);
    for (j=0; j<numpp; j++) {
        dec_ctx->otoc_mapping[j]   = j;             // original to current mapping
        dec_ctx->ctoo_mapping[j]   = j;             // current to original mapping
    }
    // Apply precoding matrix
    for (i=0; i<dec_ctx->sc->params.cnum; i++) {
        dec_ctx->JMBcoefficient[dec_ctx->sc->params.snum+dec_ctx->aoh+i][dec_ctx->sc->params.snum+i] = 1;

        NBR_node *variable_node = dec_ctx->sc->graph->l_nbrs_of_r[i]->first;        //ldpc_graph->nbrs_of_right[i];
        while (variable_node != NULL) {
            // 标记与该check packet连结的所有source packet node
            int src_pktid = variable_node->data;                        //variable_node->nb_index;
            dec_ctx->JMBcoefficient[dec_ctx->sc->params.snum+dec_ctx->aoh+i][src_pktid] = variable_node->ce;
            //dec_ctx->JMBcoefficient[NUM_SRC+OHS+i][src_pktid] = variable_node->nb_ce;
            variable_node = variable_node->next;
        }
    }


    // Step 1, translate LEVs to GEV and move them to GDM
    GF_ELEMENT *global_ces = calloc(numpp, sizeof(GF_ELEMENT));
    int p_copy = 0;                             // 拷贝到JMBcofficient的行指针
    for (i=0; i<dec_ctx->sc->params.gnum; i++) {
        matrix = dec_ctx->Matrices[i];
        for (j=0; j<gensize; j++) {
            if (matrix->coefficient[j][j] == 0)
                continue;                       // there is no local DoF here
            else {
                memset(global_ces, 0, numpp*sizeof(GF_ELEMENT));    /* Reset before reuse */
                for (k=0; k<gensize; k++) {
                    //global_ces[matrix->indices[k]] = matrix->coefficient[j][k];
                    global_ces[dec_ctx->sc->gene[i]->pktid[k]] = matrix->coefficient[j][k];
                }
                memcpy(dec_ctx->JMBcoefficient[p_copy], global_ces, numpp*sizeof(GF_ELEMENT));
                memcpy(dec_ctx->JMBmessage[p_copy], matrix->message[j], pktsize*sizeof(GF_ELEMENT));
                p_copy += 1;
            }
        }
    }
    free(global_ces);
    if (get_loglevel() == TRACE)
        printf("%d local DoFs are available, copied %d to GDM.\n", dec_ctx->local_DoF, p_copy);
    // Free up local matrices
    for (i=0; i<dec_ctx->sc->params.gnum; i++) {
        free_running_matrix(dec_ctx->Matrices[i], dec_ctx->sc->params.size_g);
        dec_ctx->Matrices[i] == NULL;
    }
    free(dec_ctx->Matrices);
    dec_ctx->Matrices = NULL;

    /* Transform GDM to upper trianguler via pivoting */
    long long ops;
    if (getenv("SNC_OA_ONEROUND") != NULL
         && atoi(getenv("SNC_OA_ONEROUND")) == 1) {
        ops = pivot_matrix_oneround(numpp+dec_ctx->aoh, numpp, pktsize, dec_ctx->JMBcoefficient, dec_ctx->JMBmessage, dec_ctx->otoc_mapping, &(dec_ctx->inactives));
    } else {
        ops = pivot_matrix_tworound(numpp+dec_ctx->aoh, numpp, pktsize, dec_ctx->JMBcoefficient, dec_ctx->JMBmessage, dec_ctx->otoc_mapping, &(dec_ctx->inactives));
    }
    dec_ctx->operations += ops;
    dec_ctx->ops2 += ops;
    /* Count available innovative rows */
    for (i=0; i<numpp; i++) {
        if (dec_ctx->JMBcoefficient[i][i] != 0)
            dec_ctx->global_DoF++;
        dec_ctx->ctoo_mapping[dec_ctx->otoc_mapping[i]] = i;
    }
}

/**
 * Save a decoding context to a file
 * Return values:
 *   On success: bytes written
 *   On error: -1
 */
long save_dec_context_OA(struct decoding_context_OA *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = OA_DECODER;
    FILE *fp;
    if ((fp = fopen(filepath, "w")) == NULL) {
        fprintf(stderr, "Cannot open %s to save decoding context\n", filepath);
        return (-1);
    }
    int i, j;
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;
    // Write snc params
    filesize += fwrite(&dec_ctx->sc->params, sizeof(struct snc_parameters), 1, fp);
    // Write decoder type
    filesize += fwrite(&(d_type), sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->aoh, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->finished, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->OA_ready, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->local_DoF, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->global_DoF, sizeof(int), 1, fp);
    // Save running matrices
    for (i=0; i<dec_ctx->sc->params.gnum; i++) {
        for (j=0; j<gensize; j++) {
            filesize += fwrite(dec_ctx->Matrices[i]->coefficient[j], sizeof(GF_ELEMENT), gensize, fp);
            filesize += fwrite(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }

    // Save GDM and its related bookkeeping information if OA ready.
    // Fixme: GDM is sparse, so we can save space by using compact representation
    // of the matrix. For example, save [row, col, value] for nonzero elements.
    if (dec_ctx->OA_ready == 1) {
        for (i=0; i<numpp + dec_ctx->aoh; i++) {
            filesize += fwrite(dec_ctx->JMBcoefficient[i], sizeof(GF_ELEMENT), numpp, fp);
            filesize += fwrite(dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
        filesize += fwrite(dec_ctx->otoc_mapping, sizeof(int), numpp, fp);
        filesize += fwrite(dec_ctx->ctoo_mapping, sizeof(int), numpp, fp);
        filesize += fwrite(&dec_ctx->inactives, sizeof(int), 1, fp);
    }
    // Save performance index
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops1, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops2, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops3, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops4, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}

struct decoding_context_OA *restore_dec_context_OA(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_parameters sp;
    fread(&sp, sizeof(struct snc_parameters), 1, fp);
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field
    int aoh;
    fread(&aoh, sizeof(int), 1, fp);
    // Create a fresh decoding context
    struct decoding_context_OA *dec_ctx = create_dec_context_OA(&sp, aoh);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_GG failed\n");
        return NULL;
    }
    // Restore decoding context from file
    int i, j;
    fread(&dec_ctx->finished, sizeof(int), 1, fp);
    fread(&dec_ctx->OA_ready, sizeof(int), 1, fp);
    fread(&dec_ctx->local_DoF, sizeof(int), 1, fp);
    fread(&dec_ctx->global_DoF, sizeof(int), 1, fp);
    // Restore running matrices
    // Note that running matrices' memory were already allocated in creating_dec_context
    for (i=0; i<dec_ctx->sc->params.gnum; i++) {
        for (j=0; j<sp.size_g; j++) {
            fread(dec_ctx->Matrices[i]->coefficient[j], sizeof(GF_ELEMENT), sp.size_g, fp);
            fread(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), sp.size_p, fp);
        }
    }
    // Restore GDM and its related information if OA_ready
    int numpp = dec_ctx->sc->params.snum + dec_ctx->sc->params.cnum;
    if (dec_ctx->OA_ready == 1) {
        dec_ctx->JMBcoefficient = calloc(numpp+aoh, sizeof(GF_ELEMENT*));
        dec_ctx->JMBmessage = calloc(numpp+aoh, sizeof(GF_ELEMENT*));
        for (i=0; i<numpp+aoh; i++) {
            dec_ctx->JMBcoefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
            fread(dec_ctx->JMBcoefficient[i], sizeof(GF_ELEMENT), numpp, fp);
            dec_ctx->JMBmessage[i] = calloc(sp.size_p, sizeof(GF_ELEMENT));
            fread(dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT), sp.size_p, fp);
        }
        dec_ctx->otoc_mapping = calloc(numpp, sizeof(int));
        fread(dec_ctx->otoc_mapping, sizeof(int), numpp, fp);
        dec_ctx->ctoo_mapping = calloc(numpp, sizeof(int));
        fread(dec_ctx->ctoo_mapping, sizeof(int), numpp, fp);
        fread(&dec_ctx->inactives, sizeof(int), 1, fp);
    }
    // Restore performance index
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops1, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops2, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops3, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops4, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}
