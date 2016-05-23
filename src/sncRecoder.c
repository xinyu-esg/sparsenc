#include <math.h>
#include "common.h"
#include "galois.h"
#include "sparsenc.h"

/* Schedule a subgeneration to recode a packet according
 * to the specified scheduling type. */
static int schedule_recode_generation(struct snc_buffer *buf, int sched_t);
static int banded_nonuniform_sched(struct snc_buffer *buf);

struct snc_buffer *snc_create_buffer(struct snc_parameters *sp, int bufsize)
{
    static char fname[] = "snc_create_decoder";
    int i;
    struct snc_buffer *buf;
    if ((buf = calloc(1, sizeof(struct snc_buffer))) == NULL) {
        fprintf(stderr, "calloc snc_buffer fail\n");
        goto Error;
    }
    //initialize recoding buffer
    buf->params = *sp;
    // determine number of generations with sp
    int num_src = ALIGN(buf->params.datasize, buf->params.size_p);
    buf->snum = num_src;
    int num_chk;
    // FIXME: such env variables should be contained in snc_parameters because
    // the settings need to be delivered to network codes in reality. For simulations
    // in which source and intermediate nodes usually are on the same host, in the same
    // process, using env variables is OK.
    char *pc = getenv("SNC_PRECISE_CHECK");
    if ( pc != NULL && atoi(pc) == 1)
        num_chk = (int) ceil(num_src * buf->params.pcrate);
    else
        num_chk = number_of_checks(num_src, buf->params.pcrate);
    buf->cnum = num_chk;
    if (buf->params.type == BAND_SNC)
        buf->gnum  = ALIGN((num_src+num_chk-buf->params.size_g), buf->params.size_b) + 1;
    else
        buf->gnum  = ALIGN( (num_src+num_chk), buf->params.size_b);

    buf->size = bufsize;
    buf->nemp = 0;
    if ((buf->gbuf = calloc(buf->gnum, sizeof(struct snc_packet **))) == NULL) {
        fprintf(stderr, "%s: calloc buf->gbuf\n", fname);
        goto Error;
    }
    for (i=0; i<buf->gnum; i++) {
        /* Initialize pointers of buffered packets of each generation as NULL */
        if ((buf->gbuf[i] = calloc(bufsize, sizeof(struct snc_packet *))) == NULL) {
            fprintf(stderr, "%s: calloc buf->gbuf[%d]\n", fname, i);
            goto Error;
        }
    }
    if ((buf->nc = calloc(buf->gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->nc\n", fname);
        goto Error;
    }
    if ((buf->pn = calloc(buf->gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->pn\n", fname);
        goto Error;
    }
    if ((buf->nsched = calloc(buf->gnum, sizeof(int))) == NULL) {
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

    struct snc_packet *pkt = snc_alloc_empty_packet(&buf->params);
    if (pkt == NULL)
        return NULL;

    pkt->gid = gid;
    GF_ELEMENT co = 0;
    for (int i=0; i<buf->nc[gid]; i++) {
        /* Coding coefficients of bufferred packets are coded together */
        if (buf->params.bnc == 1) {
            co = rand() % 2;
            if (co == 1)
                galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, ALIGN(buf->params.size_g, 8));
        } else {
            co = rand() % (1 << 8);
            galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, buf->params.size_g);
        }
        galois_multiply_add_region(pkt->syms, buf->gbuf[gid][i]->syms, co, buf->params.size_p);
    }
    return pkt;
}

int snc_recode_packet_im(struct snc_buffer *buf, struct snc_packet *pkt, int sched_t)
{
    int gid = schedule_recode_generation(buf, sched_t);
    if (gid == -1)
        return -1;
    pkt->gid = gid;
    // Clean up pkt
    if (buf->params.bnc) {
        memset(pkt->coes, 0, ALIGN(buf->params.size_g, 8)*sizeof(GF_ELEMENT));
    } else {
        memset(pkt->coes, 0, buf->params.size_g*sizeof(GF_ELEMENT));
    }
    memset(pkt->syms, 0, sizeof(GF_ELEMENT)*buf->params.size_p);
    GF_ELEMENT co = 0;
    for (int i=0; i<buf->nc[gid]; i++) {
        if (buf->params.bnc == 1) {
            co = rand() % 2;
            if (co == 1)
                galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, ALIGN(buf->params.size_g, 8));
        } else {
            co = rand() % (1 << 8);
            galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, buf->params.size_g);
        }
        galois_multiply_add_region(pkt->syms, buf->gbuf[gid][i]->syms, co, buf->params.size_p);
    }
    return 0;
}

void snc_free_buffer(struct snc_buffer *buf)
{
    if (buf == NULL)
        return;
    for (int i=0; i<buf->gnum; i++) {
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
    if (buf->nemp == 0)
        return -1;

    int gid;
    if (sched_t == TRIV_SCHED) {
        gid = rand() % buf->gnum;
        buf->nsched[gid]++;
        return gid;
    }

    if (sched_t == RAND_SCHED) {
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

    if (sched_t == NURAND_SCHED) {
        return banded_nonuniform_sched(buf);
    }
}


/*
 * Non-uniform random scheduling for banded codes
 * NOTE: scheduling of the 0-th and the (M-G)-th generation are not uniform
 * 0-th and (M-G)-th: (G+1)/2M
 * 1-th to (M-G-1)-th: 1/M
 * [G+1, 2, 2, 2,..., 2, G+1]
 * [-----{  2*(M-G-1)  }----]
 */
static int banded_nonuniform_sched(struct snc_buffer *buf)
{
	int M = buf->snum + buf->cnum;
	int G = buf->params.size_g;
	int upperb = 2*(G+1)+2*(M-G-1);

	int found = 0;
	int selected = -1;
	while (found ==0) {
		selected = (rand() % upperb) + 1;

		if (selected <= G+1) {
			selected = 0;
		} else if (selected > (G+1+2*(M-G-1))){
			selected = buf->gnum - 1;
		} else {
			int residual = selected - (G+1);
			int mod = residual / 2;
			selected = mod + 1;
		}
		if (buf->nc[selected] != 0)
			found = 1;
	}
	return selected;
}

