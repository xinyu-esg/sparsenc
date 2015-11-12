#include "galois.h"
#include "slncRecoder.h"

static int schedule_recode_generation(struct slnc_buffer *buf, int sched_t);

int slnc_create_recoding_context(struct slnc_recoding_context *rc, struct slnc_metainfo meta, int bufsize)
{
    static char fname[] = "slnc_create_recoding_context";
    int i;
    rc->meta = meta;
    //initialize recoding buffer
    rc->buf.gnum = meta.gnum;
    rc->buf.size = bufsize;
    rc->buf.nemp = 0;
    if ((rc->buf.gbuf = calloc(meta.gnum, sizeof(struct slnc_packet **))) == NULL) {
        fprintf(stderr, "%s: calloc recoding_context\n", fname);
        goto Error;
    }
    for (i=0; i<meta.gnum; i++) {
        /* Initialize pointers of buffered packets of each generation as NULL */
        if ((rc->buf.gbuf[i] = calloc(bufsize, sizeof(struct slnc_packet *))) == NULL) {
            fprintf(stderr, "%s: calloc buf.gbuf[%d]\n", fname, i);
            goto Error;
        }
    }
    if ((rc->buf.nc = calloc(meta.gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf.nc\n", fname);
        goto Error;
    }
    if ((rc->buf.pn = calloc(meta.gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf.pn\n", fname);
        goto Error;
    }
    if ((rc->buf.nsched = calloc(meta.gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf.nsched\n", fname);
        goto Error;
    }
    return (0);

Error:
    slnc_free_recoding_buffer(rc);
    return (-1);
}

/*
 * Buffer structure example (nc=1, pn=1)
 * slnc_packet     
 *      ^           NULL         NULL
 *      |            |            |
 *      |            |            |
 * gbuf[gid][0] gbuf[gid][1] gbuf[gid][2] ......... gbuf[gid][size-1] 
 *                   ^
 *                   |
 *                   |
 *                pn = 1
 */
void slnc_buffer_packet(struct slnc_recoding_context *rc, struct slnc_packet *pkt)
{
    int gid = pkt->gid;
    if (rc->buf.nc[gid] == 0) {
        /* Buffer of the generation is empty */
        rc->buf.gbuf[gid][0] = pkt;
        rc->buf.nc[gid]++;
        rc->buf.nemp++;
    } else if (rc->buf.nc[gid] == rc->buf.size) {
        /* Buffer of the generation is full, FIFO */
        slnc_free_packet(rc->buf.gbuf[gid][rc->buf.pn[gid]]);    /*discard packet previously stored in the position*/
        rc->buf.gbuf[gid][rc->buf.pn[gid]] = pkt;
    } else {
        /* Buffer is neither empty nor full */
        rc->buf.gbuf[gid][rc->buf.pn[gid]] = pkt;
        rc->buf.nc[gid]++;
    }
    rc->buf.pn[gid] = (rc->buf.pn[gid] + 1) % rc->buf.size;     /*update position for next incoming packet*/
    return;
}

struct slnc_packet *slnc_generate_recoded_packet(struct slnc_recoding_context *rc, int sched_t)
{
    int gid = schedule_recode_generation(&rc->buf, sched_t);
    if (gid == -1) 
        return NULL;

    struct slnc_packet *pkt = slnc_alloc_empty_packet(rc->meta.size_g, rc->meta.size_p);
    if (pkt == NULL)
        return NULL;

    pkt->gid = gid;
    for (int i=0; i<rc->buf.nc[gid]; i++) {
        GF_ELEMENT co = rand() % (1 << GF_POWER);
        /* Coding coefficients of bufferred packets are coded together */
        galois_multiply_add_region(pkt->coes, rc->buf.gbuf[gid][i]->coes, co, rc->meta.size_g, GF_POWER);
        galois_multiply_add_region(pkt->syms, rc->buf.gbuf[gid][i]->syms, co, rc->meta.size_p, GF_POWER);
    }
    return pkt;
}

static int schedule_recode_generation(struct slnc_buffer *buf, int sched_t)
{
    int gid;
    if (sched_t == TRIV_SCHED) {
        gid = rand() % buf->gnum;
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
        for (int j=0; j<buf->gnum; j++) {
            if (buf->nc[j] - buf->nsched[j] > max) {
                max = buf->nc[j] - buf->nsched[j];
                gid = j;
            }
        }
        buf->nsched[gid]++;
        return gid;
    }
}

void slnc_free_recoding_buffer(struct slnc_recoding_context *rc)
{
    for (int i=0; i<rc->meta.gnum; i++) {
        if (rc->buf.gbuf != NULL && rc->buf.gbuf[i] != NULL) {
            /* Free bufferred packets, if any */
            for (int j=0; j<rc->buf.size; j++)
                slnc_free_packet(rc->buf.gbuf[i][j]);	
            /* Free the pointer array */
            free(rc->buf.gbuf[i]);
        }
    }
    if (rc->buf.gbuf != NULL)
        free(rc->buf.gbuf);
    if (rc->buf.nc != NULL)
        free(rc->buf.nc);
    if (rc->buf.pn != NULL)
        free(rc->buf.pn);
    if (rc->buf.nsched != NULL)
        free(rc->buf.nsched);

    return;
}
