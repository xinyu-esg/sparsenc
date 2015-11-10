#include "gncEncoder.h"

#define TRIV_SCHED    0      // trivial scheduling
#define RAND_SCHED    1      // random scheduling a generation to generate a recoded packet
#define MLPI_SCHED    2      // Maximum Local Potential Innovativeness scheduling

/* 
 * Buffer for storing GNC packets
 *
 * Buffer size specifies how many packets are saved for 
 * each generation. "FIFO" strategy is used when buffer
 * size is reached, i.e., the oldest buffered packet is
 * discarded when a generation buffer is full while a new
 * packet belonging to the generation arrives.
 *
 * Buffer data structure
 *
 * gbuf --> gbuf[0]
 *                     coded_packet  coded_packet ...
 *          gbuf[1]          ^            ^
 *                           |            |
 *          gbuf[2] --> gbuf[2][0]   gbuf[2][1] ....
 *            .             
 *            .            
 *            .            
 */
struct gnc_buffer {
    int                    gnum;     // Number of generations
    int                    size;     // Number of bufferred packets of each generation
    int                    nemp;     // Number of non-empty buffers (i.e., have at least one packet
    // belonging to the generation is received
    struct coded_packet ***gbuf;     // Pointers to generation buffers
    int                   *nc;       // Number of currently buffered packets
    int                   *pn;       // Position where to store a new arrived packet in each generation buffer 
    int                   *nsched;   // Number of scheduled times of each generation (not used in RAND_SCHED)
};

struct gnc_recoding_context {
    struct gnc_metainfo  meta;
    struct gnc_buffer    buf;
};

int create_recoding_context(struct gnc_recoding_context *rc, struct gnc_metainfo meta, int bufsize);
void buffer_packet(struct gnc_recoding_context *rc, struct coded_packet *pkt);
struct coded_packet *generate_recoded_packet(struct gnc_recoding_context *rc, int sched_t);
void free_recoding_buffer(struct gnc_recoding_context *rc);
