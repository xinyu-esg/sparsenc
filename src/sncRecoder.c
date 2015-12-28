#include "common.h"
#include "galois.h"
#include "snc.h"

/* Schedule a subgeneration to recode a packet according
 * to the specified scheduling type. */
static int schedule_recode_generation(struct snc_buffer *buf, int sched_t);

struct snc_buffer *snc_create_buffer(struct snc_metainfo meta, int bufsize)
{
    static char fname[] = "snc_create_decoder";
    int i;
    struct snc_buffer *buf;
    if ((buf = calloc(1, sizeof(struct snc_buffer))) == NULL) {
        fprintf(stderr, "calloc snc_buffer fail\n");
        goto Error;
    }
    //initialize recoding buffer
    buf->meta = meta;
    buf->size = bufsize;
    buf->nemp = 0;
    if ((buf->gbuf = calloc(meta.gnum, sizeof(struct snc_packet **))) == NULL) {
        fprintf(stderr, "%s: calloc buf->gbuf\n", fname);
        goto Error;
    }
    for (i=0; i<meta.gnum; i++) {
        /* Initialize pointers of buffered packets of each generation as NULL */
        if ((buf->gbuf[i] = calloc(bufsize, sizeof(struct snc_packet *))) == NULL) {
            fprintf(stderr, "%s: calloc buf->gbuf[%d]\n", fname, i);
            goto Error;
        }
    }
    if ((buf->nc = calloc(meta.gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->nc\n", fname);
        goto Error;
    }
    if ((buf->pn = calloc(meta.gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->pn\n", fname);
        goto Error;
    }
    if ((buf->nsched = calloc(meta.gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->nsched\n", fname);
        goto Error;
    }
    return buf;

Error:
    snc_free_buffer(buf);
    return NULL;
}

/*
 * Buffer structure example (nc=1, pn=1)
 * snc_packet
 *      ^           NULL         NULL
 *      |            |            |
 *      |            |            |
 * gbuf[gid][0] gbuf[gid][1] gbuf[gid][2] ......... gbuf[gid][size-1]
 *                   ^
 *                   |
 *                   |
 *                pn = 1
 */
void snc_buffer_packet(struct snc_buffer *buf, struct snc_packet *pkt)
{
    int gid = pkt->gid;
    if (buf->nc[gid] == 0) {
        /* Buffer of the generation is empty */
        buf->gbuf[gid][0] = pkt;
        buf->nc[gid]++;
        buf->nemp++;
    } else if (buf->nc[gid] == buf->size) {
        /* Buffer of the generation is full, FIFO */
        snc_free_packet(buf->gbuf[gid][buf->pn[gid]]);  //discard packet previously stored in the position
        buf->gbuf[gid][buf->pn[gid]] = pkt;
    } else {
        /* Buffer is neither empty nor full */
        buf->gbuf[gid][buf->pn[gid]] = pkt;
        buf->nc[gid]++;
    }
    buf->pn[gid] = (buf->pn[gid] + 1) % buf->size;  //update position for next incoming packet
    return;
}

struct snc_packet *snc_recode_packet(struct snc_buffer *buf, int sched_t)
{
    int gid = schedule_recode_generation(buf, sched_t);
    if (gid == -1)
        return NULL;

    struct snc_packet *pkt = snc_alloc_empty_packet(&buf->meta);
    if (pkt == NULL)
        return NULL;

    pkt->gid = gid;
    for (int i=0; i<buf->nc[gid]; i++) {
        GF_ELEMENT co = rand() % (1 << GF_POWER);
        /* Coding coefficients of bufferred packets are coded together */
        galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, buf->meta.size_g, GF_POWER);
        galois_multiply_add_region(pkt->syms, buf->gbuf[gid][i]->syms, co, buf->meta.size_p, GF_POWER);
    }
    return pkt;
}

void snc_free_buffer(struct snc_buffer *buf)
{
    if (buf == NULL)
        return;
    for (int i=0; i<buf->meta.gnum; i++) {
        if (buf->gbuf != NULL && buf->gbuf[i] != NULL) {
            /* Free bufferred packets, if any */
            for (int j=0; j<buf->size; j++)
                snc_free_packet(buf->gbuf[i][j]);
            /* Free the pointer array */
            free(buf->gbuf[i]);
        }
    }
    if (buf->gbuf != NULL)
        free(buf->gbuf);
    if (buf->nc != NULL)
        free(buf->nc);
    if (buf->pn != NULL)
        free(buf->pn);
    if (buf->nsched != NULL)
        free(buf->nsched);
    free(buf);
    buf = NULL;
    return;
}

static int schedule_recode_generation(struct snc_buffer *buf, int sched_t)
{
    int gid;
    if (sched_t == TRIV_SCHED) {
        gid = rand() % buf->meta.gnum;
        buf->nsched[gid]++;
        return gid;
    }

    if (sched_t == RAND_SCHED) {
        if (buf->nemp == 0)
            return -1;
        int index = rand() % buf->nemp;
        int i = -1;
        gid = 0;
        while ( i != index) {
            if (buf->nc[gid++] != 0)
                i++;
        }
        buf->nsched[gid-1]++;
        return gid-1;
    }

    if (sched_t == MLPI_SCHED) {
        gid = 0;
        int max = buf->nc[gid] - buf->nsched[gid];
        for (int j=0; j<buf->meta.gnum; j++) {
            if (buf->nc[j] - buf->nsched[j] > max) {
                max = buf->nc[j] - buf->nsched[j];
                gid = j;
            }
        }
        buf->nsched[gid]++;
        return gid;
    }
}

