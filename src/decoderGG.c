/*--------------------------sncGGDecoder.c---------------------
 * Implementation of generation-by-generation decoding.
 *-------------------------------------------------------------*/
#include "common.h"
#include "galois.h"
#include "decoderGG.h"

// to store matrices in processing (needed by the decoder)
struct running_matrix {
    int remaining_rows;                             // record how many linearly independent encoding vectors we have
    int remaining_cols;                             // record how many source packets remain unknown
    FLAGS erased;                                   // bits indicating erased columns due to back-substitution
    GF_ELEMENT **coefficient;
    GF_ELEMENT **message;
};


static void decode_generation(struct decoding_context_GG *dec_ctx, int gid);
static void perform_iterative_decoding(struct decoding_context_GG *dec_ctx);
static void new_decoded_source_packet(struct decoding_context_GG *dec_ctx, int pkt_id);
static void new_decoded_check_packet(struct decoding_context_GG *dec_ctx, int pkt_id);
static void update_generations(struct decoding_context_GG *dec_ctx);
static long update_running_matrix(struct decoding_context_GG *dec_ctx, int gid, int sid, int index);
static int check_for_new_recoverables(struct decoding_context_GG *dec_ctx);
static int check_for_new_decodables(struct decoding_context_GG *dec_ctx);
static void mask_packet(struct decoding_context_GG *dec_ctx, GF_ELEMENT ce, int index, struct snc_packet *enc_pkt);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);

// setup decoding context:
struct decoding_context_GG *create_dec_context_GG(struct snc_parameter sp)
{
    static char fname[] = "snc_create_dec_context_GG";
    int i, j;

    struct decoding_context_GG *dec_ctx;
    if ((dec_ctx = malloc(sizeof(struct decoding_context_GG))) == NULL) {
        fprintf(stderr, "%s: malloc decoding context GG failed\n", fname);
        return NULL;
    }
    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) {
        fprintf(stderr, "%s: create decoding context failed", fname);
        goto AllocError;
    }
    dec_ctx->sc = sc;

    // memory areas needed for decoding
    dec_ctx->evolving_checks = calloc(dec_ctx->sc->meta.cnum, sizeof(GF_ELEMENT *));
    if (dec_ctx->evolving_checks == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->evolving_checks\n", fname);
        goto AllocError;
    }
    dec_ctx->check_degrees = calloc(dec_ctx->sc->meta.cnum, sizeof(int));
    if (dec_ctx->check_degrees == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->check_degrees\n", fname);
        goto AllocError;
    }
    for (i=0; i<dec_ctx->sc->meta.cnum; i++) {
        NBR_node *nb = dec_ctx->sc->graph->l_nbrs_of_r[i]->first;
        while (nb != NULL) {
            dec_ctx->check_degrees[i]++;                            /* initial check degree of each check packet */
            nb = nb->next;
        }
    }

    dec_ctx->finished  = 0;
    dec_ctx->decoded   = 0;
    dec_ctx->originals = 0;
    dec_ctx->Matrices = calloc(dec_ctx->sc->meta.gnum, sizeof(struct running_matrix*));
    if (dec_ctx->Matrices == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->Matrices\n", fname);
        goto AllocError;
    }

    for (i=0; i<dec_ctx->sc->meta.gnum; i++) {
        dec_ctx->Matrices[i] = calloc(1, sizeof(struct running_matrix));
        if (dec_ctx->Matrices[i] == NULL) {
            fprintf(stderr, "%s: malloc dec_ctx->Matrices[%d]\n", fname, i);
            goto AllocError;
        }
        dec_ctx->Matrices[i]->remaining_cols = dec_ctx->sc->meta.size_g;

        // Allocate coefficient and message matrices in running_matrix
        // coefficeint: size_g x size_g
        // message:     size_g x size_p
        // Dim-1) Pointers to each row
        dec_ctx->Matrices[i]->coefficient = calloc(dec_ctx->sc->meta.size_g, sizeof(GF_ELEMENT*));
        if (dec_ctx->Matrices[i]->coefficient == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->coefficient\n", fname, i);
            goto AllocError;
        }
        dec_ctx->Matrices[i]->message = calloc(dec_ctx->sc->meta.size_g, sizeof(GF_ELEMENT*));
        if (dec_ctx->Matrices[i]->message == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage\n", fname, i);
            goto AllocError;
        }
        // Dim-2) Elements of each row
        for (int j=0; j<dec_ctx->sc->meta.size_g; j++) {
            dec_ctx->Matrices[i]->coefficient[j] = calloc(dec_ctx->sc->meta.size_g, sizeof(GF_ELEMENT));
            if (dec_ctx->Matrices[i]->coefficient[j] == NULL) {
                fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->coefficient[%d]\n", fname, i, j);
                goto AllocError;
            }
            dec_ctx->Matrices[i]->message[j] = calloc(dec_ctx->sc->meta.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->Matrices[i]->message[j] == NULL) {
                fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage[%d]\n", fname, i, j);
                goto AllocError;
            }
        }
    }
    if ( (dec_ctx->recent = malloc(sizeof(ID_list))) == NULL ) {
        fprintf(stderr, "%s: malloc dec_ctx->recent", fname);
        goto AllocError;
    }

    dec_ctx->recent->first = dec_ctx->recent->last = NULL;
    memset(dec_ctx->grecent, -1, sizeof(int)*FB_THOLD);             /* set recent decoded generation ids to -1 */
    dec_ctx->newgpos    = 0;
    dec_ctx->grcount    = 0;
    dec_ctx->operations = 0;
    dec_ctx->overhead   = 0;
    return dec_ctx;

AllocError:
    free_dec_context_GG(dec_ctx);
    dec_ctx = NULL;
    return NULL;
}

void free_dec_context_GG(struct decoding_context_GG *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
    if (dec_ctx->sc != NULL)
        snc_free_enc_context(dec_ctx->sc);

    int i, j, k;
    if (dec_ctx->evolving_checks != NULL) {
        for (i=0; i<dec_ctx->sc->meta.cnum; i++) {
            if (dec_ctx->evolving_checks[i] != NULL)
                free(dec_ctx->evolving_checks[i]);
        }
        free(dec_ctx->evolving_checks);
    }
    if (dec_ctx->check_degrees != NULL)
        free(dec_ctx->check_degrees);
    if (dec_ctx->Matrices != NULL) {
        for (i=0; i<dec_ctx->sc->meta.gnum; i++){
            // Free each decoding matrix
            if (dec_ctx->Matrices[i] != NULL) {
                if (dec_ctx->Matrices[i]->coefficient != NULL) {
                    for (j=0; j<dec_ctx->sc->meta.size_g; j++) {
                        if (dec_ctx->Matrices[i]->coefficient[j] != NULL)
                            free(dec_ctx->Matrices[i]->coefficient[j]);
                    }
                    free(dec_ctx->Matrices[i]->coefficient);
                }
                if (dec_ctx->Matrices[i]->message != NULL) {
                    for (j=0; j<dec_ctx->sc->meta.size_g; j++) {
                        if (dec_ctx->Matrices[i]->message[j] != NULL)
                            free(dec_ctx->Matrices[i]->message[j]);
                    }
                    free(dec_ctx->Matrices[i]->message);
                }
                free(dec_ctx->Matrices[i]);
            }
        }
        free(dec_ctx->Matrices);
    }
    if (dec_ctx->recent != NULL)
        free_list(dec_ctx->recent);
    free(dec_ctx);
    dec_ctx = NULL;
    return;
}


// cache a new received packet, extract its information and try decode the class it belongs to
void process_packet_GG(struct decoding_context_GG *dec_ctx, struct snc_packet *pkt)
{
    dec_ctx->overhead += 1;

    // 1, check if the class of this packet belongs has already decoded
    int gid = pkt->gid;

    struct running_matrix *matrix = dec_ctx->Matrices[gid];                 // get the corresponding matrix in processing
    int r_rows = matrix->remaining_rows;                                    // record rows of current matrix
    int r_cols = matrix->remaining_cols;                                    // record cols of current matrix
    if (r_cols == 0)
        return;                                                             // this class has finished decoding


    // 2, extract its information, mask it against decoded packets if it's necessary
    int i = 0, j, k;
    for (j=0; j<dec_ctx->sc->meta.size_g; j++) {
        GF_ELEMENT coe;
        if (dec_ctx->sc->meta.bnc) {
            coe = get_bit_in_array(pkt->coes, j);
        } else {
            coe = pkt->coes[j];
        }
        if ( _FLAG_ON(matrix->erased, j) ) {
            //find the decoded packet, mask it with this source packet
            int src_id = dec_ctx->sc->gene[gid]->pktid[j];      // index of the corresponding source packet
            mask_packet(dec_ctx, coe, src_id, pkt);
            dec_ctx->operations += dec_ctx->sc->meta.size_p;
        } else {
            matrix->coefficient[r_rows][i] = coe;
            i++;
        }
    }
    memcpy(matrix->message[r_rows], pkt->syms, sizeof(GF_ELEMENT)*dec_ctx->sc->meta.size_p);
    matrix->remaining_rows += 1;
    snc_free_packet(pkt);
    pkt = NULL;

    //3, check if this class is full rank with at most MIN_DEG packets unknown
    //
    if (r_rows >= r_cols - 1) {
        // perform forward substitution
        long long flushing_ops = forward_substitute(r_rows+1, r_cols, dec_ctx->sc->meta.size_p, matrix->coefficient, matrix->message);
        dec_ctx->operations += flushing_ops;

        // 4, check if the new packet is full rank, if yes, decode it, if not refresh the coefficient and message matrix anyway
        int innovatives = 0;
        for (j=0; j<r_cols; j++)
            if (matrix->coefficient[j][j] != 0)
                innovatives++;

        if (innovatives != r_cols)
            matrix->remaining_rows = innovatives;
        else {
            matrix->remaining_rows = matrix->remaining_cols;
            decode_generation(dec_ctx, gid);
#if defined(GNCTRACE)
            printf("Entering perform_iterative_decoding...\n");
#endif
            perform_iterative_decoding(dec_ctx);
        }
    }
}

// decode packets of a generation via Gaussion elimination
static void decode_generation(struct decoding_context_GG *dec_ctx, int gid)
{
    static char fname[] = "decode_generation";
    //printf("entering decoding_class()...\n");
    struct running_matrix *matrix = dec_ctx->Matrices[gid];
    const int r_rows = matrix->remaining_rows;
    const int r_cols = matrix->remaining_cols;
#if defined(GNCTRACE)
    printf("%s %d: %d rows %d columns\n", fname, gid, r_rows, r_cols);
#endif

    int i, j, k;
    // this class have enough linearly independent encoding vectors, and can be decoded completely
    long long decoding_ops = back_substitute(r_rows, r_cols, dec_ctx->sc->meta.size_p, matrix->coefficient, matrix->message);
    dec_ctx->operations += decoding_ops;

    // contruct decoded packets
    int c = 0;                                              // record number of decoded pacekts
    for (i=0; i<r_rows; i++) {
        for (j=0; j<dec_ctx->sc->meta.size_g; j++) {
            // determine the index of each decoded packet
            if ( !_FLAG_ON(matrix->erased, j) ) {
                _FLAG_SET(matrix->erased, j);
                int src_id = dec_ctx->sc->gene[gid]->pktid[j];
                if (dec_ctx->sc->pp[src_id] != NULL) {
#if defined(GNCTRACE)
                    printf("%s: packet %d is already decoded.\n", fname, src_id);
#endif
                    continue;
                }
                if ( (dec_ctx->sc->pp[src_id] = calloc(dec_ctx->sc->meta.size_p, sizeof(GF_ELEMENT))) == NULL )
                    fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, src_id);
                memcpy(dec_ctx->sc->pp[src_id], matrix->message[i], sizeof(GF_ELEMENT)*dec_ctx->sc->meta.size_p);
                // Record the decoded packet as a recently decoded packet
                ID *new_id;
                if ( (new_id = malloc(sizeof(ID))) == NULL )
                    fprintf(stderr, "%s: malloc new ID\n", fname);
                new_id->data = src_id;
                new_id->next = NULL;
                append_to_list(dec_ctx->recent, new_id);
                c += 1;
                break;
            }
        }
    }
#if defined(GNCTRACE)
    printf("%s: %d packets decoded from generation %d\n", fname, c, gid);
#endif
    matrix->remaining_rows = matrix->remaining_cols = 0;
    /* Record recent decoded generations */
    dec_ctx->newgpos = dec_ctx->newgpos % FB_THOLD;
    dec_ctx->grecent[dec_ctx->newgpos] = gid;
    dec_ctx->newgpos++;
    dec_ctx->grcount++;
}

// This function performs iterative decoding on the precode and GNC code,
// based on the most recently decoded packets from a generation by decode_generation()
static void perform_iterative_decoding(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "perform_iterative_decoding";
    // Perform iterative precode decoding
    ID *new_decoded = dec_ctx->recent->first;
    while (new_decoded != NULL) {
        int new_id = new_decoded->data;
        if (new_id >= dec_ctx->sc->meta.snum)
            new_decoded_check_packet(dec_ctx, new_id);
        else
            new_decoded_source_packet(dec_ctx, new_id);
        check_for_new_recoverables(dec_ctx);
        new_decoded = new_decoded->next;
    }

    // Perform iterative generation decoding
    update_generations(dec_ctx);
    int new_decodable_gid = check_for_new_decodables(dec_ctx);
    if (new_decodable_gid != -1) {
        decode_generation(dec_ctx, new_decodable_gid);
        perform_iterative_decoding(dec_ctx);
    }
}

// Precedures to take when a source packet is decoded from a generation
static void new_decoded_source_packet(struct decoding_context_GG *dec_ctx, int pkt_id)
{
    static char fname[] = "new_decoded_source_packet";
    dec_ctx->decoded   += 1;
    dec_ctx->originals += 1;
    if (dec_ctx->originals == dec_ctx->sc->meta.snum)
        dec_ctx->finished = 1;
    /*
     * If the decoded packet comes from generations is an original source packet
     * i.e., 0<= index < NUM_SRC, back substitute into those known LDPC check packets
     */
    NBR_node *nb = dec_ctx->sc->graph->r_nbrs_of_l[pkt_id]->first;
    while (nb != NULL) {
        int check_id = nb->data;
        // If the corresponding check packet is not yet decoded, the evolving packet area
        // and the corresponding degree can be used to record the evolution of the packet.
        if (dec_ctx->evolving_checks[check_id] == NULL) {
            dec_ctx->evolving_checks[check_id] = calloc(dec_ctx->sc->meta.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->evolving_checks[check_id] == NULL)
                fprintf(stderr, "%s: calloc evolving_checks[%d]\n", fname, check_id);
        }
        // mask information bits
        galois_multiply_add_region(dec_ctx->evolving_checks[check_id], dec_ctx->sc->pp[pkt_id], nb->ce, dec_ctx->sc->meta.size_p, GF_POWER);
        dec_ctx->operations += dec_ctx->sc->meta.size_p;
        dec_ctx->check_degrees[check_id] -= 1;
        if (remove_from_list(dec_ctx->sc->graph->l_nbrs_of_r[check_id], pkt_id) == -1)
            fprintf(stderr, "%s: remove %d from l_nbrs_of_r[%d]\n", fname, pkt_id, check_id);

        nb = nb->next;
    }
}

// Procedures to take when a check packet is decoded from a generation
static void new_decoded_check_packet(struct decoding_context_GG *dec_ctx, int pkt_id)
{
    static char fname[] = "new_decoded_check_packet";
    dec_ctx->decoded += 1;

    int check_id = pkt_id - dec_ctx->sc->meta.snum;                         // ID of check packet
#if defined(GNCTRACE)
    printf("%s: degree of new decoded check %d is %d\n", \
            fname, pkt_id, dec_ctx->check_degrees[check_id]);
#endif
    if (dec_ctx->evolving_checks[check_id] == NULL) {
        // Evolving area is empty, meaning that no source neighbors of the check packet
        // has been decoded yet. So make a copy of the decoded check packet for later evolving
        dec_ctx->evolving_checks[check_id] = calloc(dec_ctx->sc->meta.size_p, sizeof(GF_ELEMENT));
        if (dec_ctx->evolving_checks[check_id] == NULL)
            fprintf(stderr, "%s: calloc evolving_checks[%d]\n", fname, check_id);
        memcpy(dec_ctx->evolving_checks[check_id], dec_ctx->sc->pp[pkt_id], sizeof(GF_ELEMENT)*dec_ctx->sc->meta.size_p);
    } else {
        // Some source neighbors have been decoded and therefore updated the evolving area,
        // we need to mask the actual check packet content against the evolving area.
        galois_multiply_add_region(dec_ctx->evolving_checks[check_id], dec_ctx->sc->pp[pkt_id], 1, dec_ctx->sc->meta.size_p, GF_POWER);
        dec_ctx->operations += dec_ctx->sc->meta.size_p;
    }
}
// This function is part of iterative precode decoding, which
// checks for new recoverable source/check packet after new packets
// are decoded from generations and processed accordingly.
static int check_for_new_recoverables(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "check_for_new_recoverables";
    int snum = dec_ctx->sc->meta.snum;
    int has_new_recoverable = -1;
    // check each check node
    for (int i=0; i<dec_ctx->sc->meta.cnum; i++) {
        if (dec_ctx->check_degrees[i] == 1
                && dec_ctx->sc->pp[i+snum] != NULL
                && !exist_in_list(dec_ctx->recent, i+snum)) {
            // The check packet is already decoded from some previous generations and its degree is
            // reduced to 1, meaning that it connects to a unrecovered source neighboer. Recover this
            // source neighbor.
            int src_id = dec_ctx->sc->graph->l_nbrs_of_r[i]->first->data;
            if (dec_ctx->sc->pp[src_id] != NULL ) {
#if defined(GNCTRACE)
                if (exist_in_list(dec_ctx->recent, src_id))
                    printf("%s: source packet %d is recoverable but is already in the recent list\n", fname, src_id);
                else
                    printf("%s: source packet %d is already decoded\n", fname, src_id);
#endif
                dec_ctx->check_degrees[i] = 0;
                continue;
            }
#if defined(GNCTRACE)
            printf("%s: source packet %d is recoverable from check %d\n", fname, src_id, i+snum);
#endif
            dec_ctx->sc->pp[src_id] = calloc(dec_ctx->sc->meta.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->sc->pp[src_id] == NULL)
                fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, src_id);
            if (dec_ctx->sc->graph->l_nbrs_of_r[i]->first->ce == 1)
                memcpy(dec_ctx->sc->pp[src_id], dec_ctx->evolving_checks[i], sizeof(GF_ELEMENT)*dec_ctx->sc->meta.size_p);
            else {
                GF_ELEMENT ce = galois_divide(1, dec_ctx->sc->graph->l_nbrs_of_r[i]->first->ce, GF_POWER);
                galois_multiply_add_region(dec_ctx->sc->pp[src_id], dec_ctx->evolving_checks[i], ce, dec_ctx->sc->meta.size_p, GF_POWER);
                dec_ctx->operations += dec_ctx->sc->meta.size_p + 1;
            }
            // Record the decoded packet as a recently decoded packet
            ID *new_id;
            if ( (new_id = malloc(sizeof(ID))) == NULL )
                fprintf(stderr, "%s: malloc new ID\n", fname);
            new_id->data = src_id;
            new_id->next = NULL;
            append_to_list(dec_ctx->recent, new_id);
            dec_ctx->check_degrees[i] = 0;
        }
        if (dec_ctx->sc->pp[i+snum] == NULL && dec_ctx->check_degrees[i] == 0) {
            // The check packet is recovered through its source neighbors
#if defined(GNCTRACE)
            printf("%s: check packet %d is recoverable\n", fname, i+snum);
#endif
            dec_ctx->sc->pp[i+snum] = calloc(dec_ctx->sc->meta.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->sc->pp[i+snum] == NULL)
                fprintf(stderr, "%s: calloc sc->pp[%d]", fname, i+snum);
            memcpy(dec_ctx->sc->pp[i+snum], dec_ctx->evolving_checks[i], sizeof(GF_ELEMENT)*dec_ctx->sc->meta.size_p);
            // Record a recently decoded packet
            ID *new_id;
            if ( (new_id = malloc(sizeof(ID))) == NULL )
                fprintf(stderr, "%s: malloc new ID\n", fname);

            new_id->data = i+snum;
            new_id->next = NULL;
            append_to_list(dec_ctx->recent, new_id);
        }

    }
    return has_new_recoverable;
}

// Update non-decoded generations with recently decoded packets
static void update_generations(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "update_generations";

    ID *precent = dec_ctx->recent->first;
    while (precent != NULL) {
        int src_id = precent->data;
        // Check all generations that contain this source packet
        for (int i=0; i<dec_ctx->sc->meta.gnum; i++) {
            if (dec_ctx->Matrices[i]->remaining_cols == 0)
                continue;
            int position = has_item(dec_ctx->sc->gene[i]->pktid, src_id, dec_ctx->sc->meta.size_g);
            //if (position != -1 && dec_ctx->Matrices[i]->erased[position] == 0) {
            if (position != -1 && _FLAG_OFF(dec_ctx->Matrices[i]->erased, position)) {
                // The recently decoded packet is not-yet decoded in the generation, mask it
                long ops = update_running_matrix(dec_ctx, i, src_id, position);
                dec_ctx->operations += ops;
            }
        }
        precent = precent->next;
    }
    // Clean up recently decoded packet ID list
    clear_list(dec_ctx->recent);
}

/********************************************************
 * updating running matrix with a decoded packet
 *   gid - generation id of the matrix to be updated;
 *   sid - packet id against which the matrix is updated;
 *   index - column index the packet belongs to
 *********************************************************/
static long update_running_matrix(struct decoding_context_GG *dec_ctx, int gid, int sid, int index)
{
    static char fname[] = "update_running_matrix";
    //printf("entering update_running_matrix()...\n");
    int i, j, k;
    long operations = 0;
    struct running_matrix *matrix = dec_ctx->Matrices[gid];
    int r_rows = matrix->remaining_rows;
    int r_cols = matrix->remaining_cols;

    // 1, update coefficient matrix: erase one column
    // index - the column position where the dec_pkt is in this matrix
    //
    // Calculate how many packet remains unknown in front of the position index of this class
    // move the erased column to the last column, move remaining columns behind forward
    int count = 0;
    for (i=0; i<index; i++) {
        if ( _FLAG_OFF(matrix->erased, i) )
            count += 1;
    }
    // Now the column of index "count" is the column corresponding to the decoded packet of ID "index"
    for (j=count+1; j<r_cols; j++) {
        for (k=0; k<r_rows; k++) {
            GF_ELEMENT tmp = matrix->coefficient[k][j-1];
            matrix->coefficient[k][j-1] = matrix->coefficient[k][j];
            matrix->coefficient[k][j] = tmp;
        }
    }
    // 2, update the message matrix
    for (i=0; i<r_rows; i++) {
        GF_ELEMENT ce = matrix->coefficient[i][r_cols-1];       // the encoding coefficient of the erasuing packet
        if (ce == 0)
            continue;                                       // NOTE: the matrix could high likely be sparse
        galois_multiply_add_region(&matrix->message[i][0], dec_ctx->sc->pp[sid], ce, dec_ctx->sc->meta.size_p, GF_POWER);
        operations += dec_ctx->sc->meta.size_p;
    }
    _FLAG_SET(matrix->erased, index);       // mark the column as erased
    matrix->remaining_cols -= 1;
    return operations;
}

// Check if there is new decodable generations
static int check_for_new_decodables(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "check_for_new_decodables";
    int i, j, k;
    for (i=0; i<dec_ctx->sc->meta.gnum; i++) {
        struct running_matrix *matrix = dec_ctx->Matrices[i];
        if ( matrix->remaining_cols != 0
                && matrix->remaining_rows >= matrix->remaining_cols ) {
            int r_rows = matrix->remaining_rows;
            int r_cols = matrix->remaining_cols;
            // perform forward substitution
            int flushing_ops = forward_substitute(r_rows, r_cols, dec_ctx->sc->meta.size_p, matrix->coefficient, matrix->message);
            dec_ctx->operations += flushing_ops;

            // check if the new packet is full rank, if yes, decode it, if not refresh the coefficient and message matrix anyway
            int full_rank = 1;
            int innovatives = 0;
            for (j=0; j<r_cols; j++) {
                if (matrix->coefficient[j][j] == 0)
                    full_rank = 0;
                else
                    innovatives += 1;
            }

            if (full_rank == 0)
                matrix->remaining_rows = innovatives;
            else {
                matrix->remaining_rows = matrix->remaining_cols;
#if defined(GNCTRACE)
                printf("%s: generation %d is decodable\n",fname, i);
#endif
                return i;
            }
        }
    }
    return -1;
}


// mask the encoded packet with the decoded packet list
// ce    - coefficient we use in masking
// index - index of the decoded source packet we want to mask against ENC_pkt
static void mask_packet(struct decoding_context_GG *dec_ctx, GF_ELEMENT ce, int index, struct snc_packet *enc_pkt)
{
    static char fname[] = "mask_packet";
    galois_multiply_add_region(enc_pkt->syms, dec_ctx->sc->pp[index], ce, dec_ctx->sc->meta.size_p, GF_POWER);
    dec_ctx->operations += dec_ctx->sc->meta.size_p;
    return;
}

/**
 * Save a decoding context to a file
 * Return values:
 *   On success: bytes written
 *   On error: -1
 */
long save_dec_context_GG(struct decoding_context_GG *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = GG_DECODER;
    FILE *fp;
    if ((fp = fopen(filepath, "w")) == NULL) {
        fprintf(stderr, "Cannot open %s to save decoding context\n", filepath);
        return (-1);
    }
    int i, j, k;
    int gensize = dec_ctx->sc->meta.size_g;
    int pktsize = dec_ctx->sc->meta.size_p;
    int numpp   = dec_ctx->sc->meta.snum + dec_ctx->sc->meta.cnum;
    // Write snc metainfo
    filesize += fwrite(&dec_ctx->sc->meta, sizeof(struct snc_metainfo), 1, fp);
    // Write decoder type
    filesize += fwrite(&(d_type), sizeof(int), 1, fp);
    // Save already decoded packets in dec_ctx->sc->pp
    filesize += fwrite(&dec_ctx->decoded, sizeof(int), 1, fp);
    for (i=0; i<numpp; i++) {
        if (dec_ctx->sc->pp[i] != NULL) {
            filesize += fwrite(&i, sizeof(int), 1, fp);  // pktid
            filesize += fwrite(dec_ctx->sc->pp[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    // Save evolving check packets
    int count = 0;
    for (i=0; i<dec_ctx->sc->meta.cnum; i++) {
        if (dec_ctx->evolving_checks[i] != NULL)
            count++;
    }
    filesize += fwrite(&count, sizeof(int), 1, fp);  // evolving packets non-NULL count
    for (i=0; i<dec_ctx->sc->meta.cnum; i++) {
        if (dec_ctx->evolving_checks[i] != NULL) {
            filesize += fwrite(&i, sizeof(int), 1, fp);  // check id
            filesize += fwrite(dec_ctx->evolving_checks[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    // Save check degrees
    filesize += fwrite(dec_ctx->check_degrees, sizeof(int), dec_ctx->sc->meta.cnum, fp);
    filesize += fwrite(&dec_ctx->finished, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->decoded, sizeof(int), 1, fp);  // Yes, I know. I saved this value twice!
    filesize += fwrite(&dec_ctx->originals, sizeof(int), 1, fp);
    // Save running matrices
    for (i=0; i<dec_ctx->sc->meta.gnum; i++) {
        filesize += fwrite(&dec_ctx->Matrices[i]->remaining_rows, sizeof(int), 1, fp);
        filesize += fwrite(&dec_ctx->Matrices[i]->remaining_cols, sizeof(int), 1, fp);
        filesize += fwrite(&dec_ctx->Matrices[i]->erased, sizeof(FLAGS), 1, fp);
        //for (j=0; j<gensize; j++) {
        for (j=0; j<dec_ctx->Matrices[i]->remaining_rows; j++) {
            filesize += fwrite(dec_ctx->Matrices[i]->coefficient[j], sizeof(GF_ELEMENT), gensize, fp);
            filesize += fwrite(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    // Save recent ID_list
    count = 0;
    ID *id = dec_ctx->recent->first;
    while (id != NULL) {
        count++;
        id = id->next;
    }
    filesize += fwrite(&count, sizeof(int), 1, fp);  // Number of recent IDs in the list
    id = dec_ctx->recent->first;
    while (count > 0) {
        fwrite(&id->data, sizeof(int), 1, fp);
        count--;
        id = id->next;
    }
    // Save performance index
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}

struct decoding_context_GG *restore_dec_context_GG(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_metainfo meta;
    fread(&meta, sizeof(struct snc_metainfo), 1, fp);
    struct snc_parameter sp;
    sp.datasize = meta.datasize;
    sp.pcrate = meta.pcrate;
    sp.size_b = meta.size_b;
    sp.size_g = meta.size_g;
    sp.size_p = meta.size_p;
    sp.type = meta.type;
    sp.bpc = meta.bpc;
    sp.bnc = meta.bnc;
    // Create a fresh decoding context
    struct decoding_context_GG *dec_ctx = create_dec_context_GG(sp);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_GG failed\n");
        return NULL;
    }
    // Restore decoding context from file
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field
	// Restore already decoded packets
	int i, j, k;
	fread(&dec_ctx->decoded, sizeof(int), 1, fp);
	for (i=0; i<dec_ctx->decoded; i++) {
		int pktid;
		fread(&pktid, sizeof(int), 1, fp);
		dec_ctx->sc->pp[pktid] = calloc(meta.size_p, sizeof(GF_ELEMENT));
		fread(dec_ctx->sc->pp[pktid], sizeof(GF_ELEMENT), meta.size_p, fp);
	}
	// Restore evolving packets
	int count;
	fread(&count, sizeof(int), 1, fp);
	for (i=0; i<count; i++) {
		int evoid;
		fread(&evoid, sizeof(int), 1, fp);
		dec_ctx->evolving_checks[evoid] = calloc(meta.size_p, sizeof(GF_ELEMENT));
		fread(dec_ctx->evolving_checks[evoid], sizeof(GF_ELEMENT), meta.size_p, fp);
	}
	// Restore check degrees
	fread(dec_ctx->check_degrees, sizeof(int), dec_ctx->sc->meta.cnum, fp);
	fread(&dec_ctx->finished, sizeof(int), 1, fp);
	fread(&dec_ctx->decoded, sizeof(int), 1, fp);
	fread(&dec_ctx->originals, sizeof(int), 1, fp);
	// Restore running matrices
	// Note that running matrices' memory were already allocated in creating_dec_context
	for (i=0; i<dec_ctx->sc->meta.gnum; i++) {
		fread(&dec_ctx->Matrices[i]->remaining_rows, sizeof(int), 1, fp);
		fread(&dec_ctx->Matrices[i]->remaining_cols, sizeof(int), 1, fp);
		fread(&dec_ctx->Matrices[i]->erased, sizeof(FLAGS), 1, fp);
		//for (j=0; j<meta.size_g; j++) {
		for (j=0; j<dec_ctx->Matrices[i]->remaining_rows; j++) {
			fread(dec_ctx->Matrices[i]->coefficient[j], sizeof(GF_ELEMENT), meta.size_g, fp);
			fread(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), meta.size_p, fp);
		}
	}
	// Restore recent ID_list
	fread(&count, sizeof(int), 1, fp);
	ID *new_id;
	for (i=0; i<count; i++) {
		if ( (new_id = malloc(sizeof(ID))) == NULL )
			fprintf(stderr, "malloc new ID failed\n");
		fread(&new_id->data, sizeof(int), 1, fp);
		new_id->next = NULL;
		append_to_list(dec_ctx->recent, new_id);
	}
	// Restore performance index
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}
