/**************************************************************
 * sncEncoder.c
 *
 * Functions for SNC encoding. Coded packets can be generated
 * from memory buffer or files.
 **************************************************************/
#include <math.h>
#include "common.h"
#include "galois.h"
#include "sparsenc.h"

static int create_context_from_meta(struct snc_context *sc);
static int verify_code_parameter(struct snc_metainfo *meta);
static void perform_precoding(struct snc_context *sc);
static int group_packets_rand(struct snc_context *sc);
static int group_packets_band(struct snc_context *sc);
static int group_packets_windwrap(struct snc_context *sc);
static void encode_packet(struct snc_context *sc, int gid, struct snc_packet *pkt);
static int schedule_generation(struct snc_context *sc);
/*
 * Create a GNC context containing meta information about the data to be encoded.
 *   buf      - Buffer containing bytes of data to be encoded
 *   sc       - pointer to the snc_context where the context will be stored
 *   snc_context is defined in src/common.h
 *
 * Return
 *   0  - Create successfully
 *   -1 - Create failed
 */
struct snc_context *snc_create_enc_context(unsigned char *buf, struct snc_parameter sp)
{
    static char fname[] = "snc_create_enc_context";
    // Allocate file_context
    struct snc_context *sc;
    if ( (sc = calloc(1, sizeof(struct snc_context))) == NULL ) {
        fprintf(stderr, "%s: calloc file_context\n", fname);
        return NULL;
    }
    sc->meta.datasize = sp.datasize;
    sc->meta.pcrate   = sp.pcrate;
    sc->meta.size_b   = sp.size_b;
    sc->meta.size_g   = sp.size_g;
    sc->meta.size_p   = sp.size_p;
    sc->meta.type     = sp.type;
    sc->meta.bpc      = sp.bpc;
    sc->meta.bnc      = sp.bnc;
    sc->meta.sys      = sp.sys;
    // Determine packet and generation numbers
    int num_src = ALIGN(sc->meta.datasize, sc->meta.size_p);
    int num_chk = number_of_checks(num_src, sc->meta.pcrate);
    sc->meta.snum  = num_src;  // Number of source packets
    sc->meta.cnum  = num_chk;  // Number of check packets
    if (sc->meta.type == BAND_SNC) {
        sc->meta.gnum  = ALIGN((num_src+num_chk-sc->meta.size_g), sc->meta.size_b) + 1;
    } else {
        sc->meta.gnum  = ALIGN( (num_src+num_chk), sc->meta.size_b);
    }
    /*
     * Verify code parameter
     */
    if (verify_code_parameter(&(sc->meta)) != 0) {
        fprintf(stderr, "%s: code parameter is invalid.\n", fname);
        snc_free_enc_context(sc);
        return NULL;
    }
    /*
     * Create generations, bipartite graph
     */
    if (create_context_from_meta(sc) != 0) {
        fprintf(stderr, "%s: create_context_from_meta\n", fname);
        snc_free_enc_context(sc);
        return NULL;
    }

    // Allocating pointers to data
    if ((sc->pp = calloc(sc->meta.snum+sc->meta.cnum, sizeof(GF_ELEMENT*))) == NULL) {
        fprintf(stderr, "%s: calloc sc->pp\n", fname);
        snc_free_enc_context(sc);
        return NULL;
    }

    constructField(GF_POWER);   // Construct Galois Field
    if (buf != NULL) {
        int alread = 0;
        int i;
        // Load source packets
        for (i=0; i<sc->meta.snum; i++) {
            sc->pp[i] = calloc(sc->meta.size_p, sizeof(GF_ELEMENT));
            int toread = (alread+sc->meta.size_p) <= sc->meta.datasize ? sc->meta.size_p : sc->meta.datasize-alread;
            memcpy(sc->pp[i], buf+alread, toread*sizeof(GF_ELEMENT));
            alread += toread;
        }
        // Allocate parity-check packet space
        for (i=0; i<sc->meta.cnum; i++)
            sc->pp[sc->meta.snum+i] = calloc(sc->meta.size_p, sizeof(GF_ELEMENT));
        perform_precoding(sc);
    }

    return sc;
}

inline struct snc_metainfo *snc_get_metainfo(struct snc_context *sc)
{
    if (sc == NULL) {
        return NULL;
    } else {
        return &(sc->meta);
    }
}

/*
 * Load data from file  into a snc_context.
 *   start - starting point to read file
 * It is caller's responsibility to ensure:
 *   start + sc->meta.datasize <= fileSize
 */
int snc_load_file_to_context(const char *filepath, long start, struct snc_context *sc)
{
    static char fname[] = "snc_load_file_to_context";

    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL)
        return (-1);
    fseek(fp, 0, SEEK_END);
    if ((ftell(fp) - start) < sc->meta.datasize) {
        return (-1);
    }
    fseek(fp, start, SEEK_SET);  // seek to position start
    int alread = 0;
    int i;
    for (i=0; i<sc->meta.snum; i++) {
        sc->pp[i] = calloc(sc->meta.size_p, sizeof(GF_ELEMENT));
        int toread = (alread+sc->meta.size_p) <= sc->meta.datasize ? sc->meta.size_p : sc->meta.datasize-alread;
        if (fread(sc->pp[i], sizeof(GF_ELEMENT), toread, fp) != toread) {
            fprintf(stderr, "%s: fread sc->pp[%d]\n", fname, i);
            return (-1);
        }
        alread += toread;
    }
    fclose(fp);
    // Allocate parity-check packet space
    for (i=0; i<sc->meta.cnum; i++)
        sc->pp[sc->meta.snum+i] = calloc(sc->meta.size_p, sizeof(GF_ELEMENT));
    perform_precoding(sc);
    return (0);
}

static int verify_code_parameter(struct snc_metainfo *meta)
{
    if (meta->size_b > meta->size_g) {
        fprintf(stderr, "code spmeter error: size_b > size_g\n");
        return(-1);
    }
    if (meta->size_g*meta->size_p > meta->datasize) {
        fprintf(stderr, "code spmeter error: size_b X size_p > datasize\n");
        return(-1);
    }
    return(0);
}

/*
 * Create snc context using metadata
 */
static int create_context_from_meta(struct snc_context *sc)
{
    static char fname[] = "snc_create_enc_context_meta";
    // Inintialize generation structures
    sc->gene  = calloc(sc->meta.gnum, sizeof(struct subgeneration*));
    if ( sc->gene == NULL ) {
        fprintf(stderr, "%s: malloc sc->gene\n", fname);
        return(-1);
    }
    for (int j=0; j<sc->meta.gnum; j++) {
        sc->gene[j] = malloc(sizeof(struct subgeneration));
        if ( sc->gene[j] == NULL ) {
            fprintf(stderr, "%s: malloc sc->gene[%d]\n", fname, j);
            return(-1);
        }
        sc->gene[j]->gid = -1;
        sc->gene[j]->pktid = malloc(sizeof(int)*sc->meta.size_g);       // Use malloc because pktid needs to be initialized as -1's later
        if ( sc->gene[j]->pktid == NULL ) {
            fprintf(stderr, "%s: malloc sc->gene[%d]->pktid\n", fname, j);
            return(-1);
        }
        memset(sc->gene[j]->pktid, -1, sizeof(int)*sc->meta.size_g);
    }
    sc->nccount = calloc(sc->meta.gnum, sizeof(int));
    if (sc->nccount == NULL) {
        fprintf(stderr, "%s: calloc sc->nccount\n", fname);
        return (-1);
    }

    int coverage;
    if (sc->meta.type == RAND_SNC) {
        coverage = group_packets_rand(sc);
    } else if (sc->meta.type == BAND_SNC) {
        coverage = group_packets_band(sc);
    } else if (sc->meta.type == WINDWRAP_SNC) {
        coverage = group_packets_windwrap(sc);
    }
#if defined(GNCTRACE)
    printf("Data Size: %ld\t Source Packets: %d\t Check Packets: %d\t Generations: %d\t Coverage: %d\n",sc->meta.datasize, sc->meta.snum, sc->meta.cnum, sc->meta.gnum, coverage);
#endif
    // Creating bipartite graph of the precode
    if (sc->meta.cnum != 0) {
        if ( (sc->graph = malloc(sizeof(BP_graph))) == NULL ) {
            fprintf(stderr, "%s: malloc BP_graph\n", fname);
            return (-1);
        }
        sc->graph->binaryce = sc->meta.bpc;     // If precode in GF(2), edges use 1 as coefficient
        if (create_bipartite_graph(sc->graph, sc->meta.snum, sc->meta.cnum) < 0)
            return (-1);
    }
    return(0);
}

void snc_free_enc_context(struct snc_context *sc)
{
    if (sc == NULL)
        return;
    int i;
    if (sc->pp != NULL) {
        for (i=sc->meta.snum+sc->meta.cnum-1; i>=0; i--) {
            if (sc->pp[i] != NULL) {
                free(sc->pp[i]);
                sc->pp[i] = NULL;
            }
        }
        free(sc->pp);
    }
    if (sc->gene != NULL) {
        for (i=sc->meta.gnum-1; i>=0; i--) {
            free(sc->gene[i]->pktid);  // free packet IDs
            free(sc->gene[i]);         // free generation itself
            sc->gene[i] = NULL;
        }
        free(sc->gene);
    }
    if (sc->graph != NULL)
        free_bipartite_graph(sc->graph);
    if (sc->nccount != NULL)
        free(sc->nccount);
    free(sc);
    sc = NULL;
    return;
}

unsigned char *snc_recover_data(struct snc_context *sc)
{
    static char fname[] = "snc_recover_data";
    long datasize = sc->meta.datasize;
    long alwrote = 0;
    long towrite = datasize;

    unsigned char *data;
    if ( (data = malloc(datasize)) == NULL) {
        fprintf(stderr, "%s: malloc(datasize) failed.\n", fname);
        return NULL;
    }
    memset(data, 0, sizeof(unsigned char)*datasize);
    int pc = 0;
    while (alwrote < datasize) {
        towrite = ((alwrote + sc->meta.size_p) <= datasize) ? sc->meta.size_p : datasize - alwrote;
        memcpy(data+alwrote, sc->pp[pc++], sizeof(GF_ELEMENT)*towrite);
        alwrote += towrite;
    }
    return data;
}

// Wrapper of free()
void snc_free_recovered(unsigned char *data)
{
    if (data != NULL)
        free(data);
    return;
}

/**
 * Recover data to file.
 * fp has to be opened in 'a' (append) mode
 **/
long snc_recover_to_file(const char *filepath, struct snc_context *sc)
{
    static char fname[] = "snc_recover_to_file";
    long datasize = sc->meta.datasize;
    long alwrote = 0;
    long towrite = datasize;

#if defined(GNCTRACE)
    printf("Writing to decoded file.\n");
#endif

    FILE *fp;
    if ((fp = fopen(filepath, "a")) == NULL)
        return (-1);

    int pc = 0;
    while (alwrote < datasize) {
        towrite = ((alwrote + sc->meta.size_p) <= datasize) ? sc->meta.size_p : datasize - alwrote;
        if (fwrite(sc->pp[pc], sizeof(GF_ELEMENT), towrite, fp) != towrite)
            fprintf(stderr, "%s: fwrite sc->pp[%d]\n", fname, pc);
        pc++;
        alwrote += towrite;
    }
    fclose(fp);
    return alwrote;
}

// perform systematic LDPC precoding against SRC pkt list and results in a LDPC pkt list
static void perform_precoding(struct snc_context *sc)
{
    static char fname[] = "perform_precoding";

    int i, j;
    for (i=0; i<sc->meta.cnum; i++) {
        // Encoding check packet according to the LDPC graph
        NBR_node *nb = sc->graph->l_nbrs_of_r[i]->first;
        while(nb != NULL) {
            int sid = nb->data;  // index of source packet
            // XOR information content
            galois_multiply_add_region(sc->pp[i+sc->meta.snum], sc->pp[sid], nb->ce, sc->meta.size_p, GF_POWER);
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
static int group_packets_rand(struct snc_context *sc)
{
    int num_p = sc->meta.snum + sc->meta.cnum;
    int num_g = sc->meta.gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int rotate = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        // split packets into disjoint groups
        for (j=0; j<sc->meta.size_b; j++) {
            index = (i * sc->meta.size_b + j) % num_p;  // source packet index

            while (has_item(sc->gene[i]->pktid, index, j) != -1)
                index++;
            sc->gene[i]->pktid[j] = index;
            selected[index] += 1;
        }

        // fill in the rest of the generation with packets from other generations
        for (j=sc->meta.size_b; j<sc->meta.size_g; j++) {
            int magicX = sc->meta.size_b + num_g - sc->meta.size_g;
            magicX = 7 <= magicX ? 7 : magicX;      // Make sure start won't be negative
            int tmp = i - (j - sc->meta.size_b + magicX);
            int start = tmp >= 0 ? tmp : tmp+num_g;
            if (start == i)
                start++;
            index = (start * sc->meta.size_b + (j - sc->meta.size_b + rotate) % (sc->meta.size_g)) % num_p;
            while (has_item(sc->gene[i]->pktid, index, j) != -1) {
                index++;
                index = index % num_p;
            }
            sc->gene[i]->pktid[j] = index;
            selected[index] += 1;
        }
        rotate = (rotate + 7) % (sc->meta.size_g);
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
static int group_packets_band(struct snc_context *sc)
{
    int num_p = sc->meta.snum + sc->meta.cnum;
    int num_g = sc->meta.gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int leading_pivot = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        leading_pivot = i * sc->meta.size_b;
        if (leading_pivot > num_p - sc->meta.size_g) {
#if defined(GNCTRACE)
            printf("Band lead of gid: %d is modified\n", i);
#endif
            leading_pivot = num_p - sc->meta.size_g;
        }
        for (j=0; j<sc->meta.size_g; j++) {
            index = leading_pivot + j;
            selected[index] += 1;
            sc->gene[i]->pktid[j] = index;
        }
    }
    int coverage = 0;
    for (i=0; i<num_p; i++)
        coverage += selected[i];

    free(selected);
    return coverage;
}

/*
 * Group packets to generations that overlap consecutively. Wrap around if needed.
 */
static int group_packets_windwrap(struct snc_context *sc)
{
    int num_p = sc->meta.snum + sc->meta.cnum;
    int num_g = sc->meta.gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int leading_pivot = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        leading_pivot = i * sc->meta.size_b;
        for (j=0; j<sc->meta.size_g; j++) {
            index = (leading_pivot + j) % num_p;
            selected[index] += 1;
            sc->gene[i]->pktid[j] = index;
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
struct snc_packet *snc_alloc_empty_packet(struct snc_metainfo *meta)
{
    struct snc_packet *pkt = calloc(1, sizeof(struct snc_packet));
    if (pkt == NULL)
        return NULL;
    if (meta->bnc) {
        pkt->coes = calloc(ALIGN(meta->size_g, 8), sizeof(GF_ELEMENT));    // Each unsigned char contains 8 bits
    } else {
        pkt->coes = calloc(meta->size_g, sizeof(GF_ELEMENT));
    }
    if (pkt->coes == NULL)
        goto AllocErr;
    pkt->syms = calloc(meta->size_p, sizeof(GF_ELEMENT));
    if (pkt->syms == NULL)
        goto AllocErr;

    return pkt;

AllocErr:
    snc_free_packet(pkt);
    return NULL;
}

/* Generate a GNC coded packet. Memory is allocated in the function. */
struct snc_packet *snc_generate_packet(struct snc_context *sc)
{
    struct snc_packet *pkt = snc_alloc_empty_packet(&sc->meta);
    int gid = schedule_generation(sc);
    encode_packet(sc, gid, pkt);
    return pkt;
}

/*
 * Generate a GNC coded packet in a given memory area.
 * It is the caller's responsibity to allocate memory properly.
 */
int snc_generate_packet_im(struct snc_context *sc, struct snc_packet *pkt)
{
    if (pkt == NULL || pkt->coes == NULL || pkt->syms == NULL)
        return -1;
    if (sc->meta.bnc) {
        memset(pkt->coes, 0, ALIGN(sc->meta.size_g, 8)*sizeof(GF_ELEMENT));
    } else {
        memset(pkt->coes, 0, sc->meta.size_g*sizeof(GF_ELEMENT));
    }
    memset(pkt->syms, 0, sc->meta.size_p*sizeof(GF_ELEMENT));
    int gid = schedule_generation(sc);
    encode_packet(sc, gid, pkt);
    return (0);
}

void snc_free_packet(struct snc_packet *pkt)
{
    if (pkt == NULL)
        return;
    if (pkt->coes != NULL)
        free(pkt->coes);
    if (pkt->syms != NULL)
        free(pkt->syms);
    free(pkt);
}


static void encode_packet(struct snc_context *sc, int gid, struct snc_packet *pkt)
{
    pkt->gid = gid;
    int pktid;
    if (sc->meta.sys == 1 && sc->nccount[gid] < sc->meta.size_b) {
        // Send an uncoded packet
        if (sc->meta.bnc == 1) {
            set_bit_in_array(pkt->coes, sc->nccount[gid]);
        } else {
            pkt->coes[sc->nccount[gid]] = 1;
        }
        pktid = sc->gene[gid]->pktid[sc->nccount[gid]];
        memcpy(pkt->syms, sc->pp[pktid], sc->meta.size_p*sizeof(GF_ELEMENT));
        sc->nccount[gid] += 1;
        return;
    }
    int i;
    GF_ELEMENT co;
    for (i=0; i<sc->meta.size_g; i++) {
        pktid = sc->gene[gid]->pktid[i];  // The i-th packet of the gid-th generation
        if (sc->meta.bnc) {
            co = (GF_ELEMENT) rand() % 2;                   // Binary network code
            if (co == 1)
                set_bit_in_array(pkt->coes, i);             // Set the corresponding coefficient as 1
        } else {
            co = (GF_ELEMENT) rand() % (1 << GF_POWER);     // Randomly generated coding coefficient
            pkt->coes[i] = co;
        }
        galois_multiply_add_region(pkt->syms, sc->pp[pktid], co, sc->meta.size_p, GF_POWER);
    }
    sc->nccount[gid] += 1;
    return;
}

static int schedule_generation(struct snc_context *sc)
{
    if (sc->meta.gnum == 1)
        return 0;
    int gid = rand() % (sc->meta.gnum);
    return gid;
}

/*
 * Print code summary
 * If called by decoders, it prints overhead and operations as well.
 */
void print_code_summary(struct snc_context *sc, int overhead, long long operations)
{
    char typestr[20];
    switch(sc->meta.type) {
        case RAND_SNC:
            strcpy(typestr, "RAND");
            break;
        case BAND_SNC:
            strcpy(typestr, "BAND");
            break;
        case WINDWRAP_SNC:
            strcpy(typestr, "WINDWRAP");
            break;
        default:
            strcpy(typestr, "UNKNOWN");
    }
    char typestr2[20];
    if (sc->meta.bpc) {
        strcpy(typestr2, "BinaryPrecode");
    } else {
        strcpy(typestr2, "NonBinaryPrecode");
    }
    char typestr3[20];
    if (sc->meta.bnc) {
        strcpy(typestr3, "BinaryNC");
    } else {
        strcpy(typestr3, "NonBinaryNC");
    }
    char typestr4[20];
    if (sc->meta.sys) {
        strcpy(typestr4, "Systematic");
    } else {
        strcpy(typestr4, "NonSystematic");
    }
    printf("datasize: %d ", sc->meta.datasize);
    printf("precode: %.3f ", sc->meta.pcrate);
    printf("size_b: %d ", sc->meta.size_b);
    printf("size_g: %d ", sc->meta.size_g);
    printf("size_p: %d ", sc->meta.size_p);
    printf("type: [%s::%s::%s::%s] ", typestr, typestr2, typestr3, typestr4);
    printf("snum: %d ", sc->meta.snum);
    printf("cnum: %d ", sc->meta.cnum);
    printf("gnum: %d ", sc->meta.gnum);
    if (operations != 0) {
        printf("overhead: %.3f ", (double) overhead/sc->meta.snum);
        printf("computation: %f\n", (double) operations/sc->meta.snum/sc->meta.size_p);
    } else {
        printf("\n");
    }
}

