#ifndef SNC_RECODER
#define SNC_RECODER
#include "sncEncoder.h"

#define TRIV_SCHED    0      // trivial scheduling
#define RAND_SCHED    1      // random scheduling a subgeneration to generate a recoded packet
#define MLPI_SCHED    2      // Maximum Local Potential Innovativeness scheduling
struct snc_metainfo;
/* 
 * Buffer for storing GNC packets
 *
 * Buffer size specifies how many packets are saved for 
 * each subgeneration. "FIFO" strategy is used when buffer
 * size is reached, i.e., the oldest buffered packet is
 * discarded when a subgeneration buffer is full while a new
 * packet belonging to the subgeneration arrives.
 *
 * Buffer data structure
 *
 * gbuf --> gbuf[0]
 *                     snc_packet  snc_packet ...
 *          gbuf[1]          ^            ^
 *                           |            |
 *          gbuf[2] --> gbuf[2][0]   gbuf[2][1] ....
 *            .             
 *            .            
 *            .            
 */
struct snc_buffer {
    int                    gnum;     // Number of subgenerations
    int                    size;     // Number of bufferred packets of each subgeneration
    int                    nemp;     // Number of non-empty buffers (i.e., have at least one packet
    // belonging to the subgeneration is received
    struct snc_packet   ***gbuf;     // Pointers to subgeneration buffers
    int                   *nc;       // Number of currently buffered packets
    int                   *pn;       // Position where to store a new arrived packet in each subgeneration buffer 
    int                   *nsched;   // Number of scheduled times of each subgeneration (not used in RAND_SCHED)
};

struct snc_recoding_context {
    struct snc_metainfo  meta;
    struct snc_buffer    buf;
};

int snc_create_recoding_context(struct snc_recoding_context *rc, struct snc_metainfo meta, int bufsize);
void snc_buffer_packet(struct snc_recoding_context *rc, struct snc_packet *pkt);
struct snc_packet *snc_generate_recoded_packet(struct snc_recoding_context *rc, int sched_t);
void snc_free_recoding_buffer(struct snc_recoding_context *rc);
#endif
