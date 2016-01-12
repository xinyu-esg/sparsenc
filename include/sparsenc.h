#ifndef SNC_H
#define SNC_H
#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
/*
 * Type of SNC codes
 * RAND     - Packets are pesudo-randomly grouped
 * BAND     - Packets are grouped to be consecutively overlapped
 * WINDWRAP - Similar as BAND, but has wrap around in encoding vectors
 */
#define RAND_SNC        0
#define BAND_SNC        1
#define WINDWRAP_SNC    2

/*
 * Type of SNC decoders
 * GG  - (Sub)Generation-by-(sub)generation decoder
 * OA  - Overlap aware decoder
 * BD  - Band decoder with pivoting
 * CBD - Compact band decoder with compact decoding matrix representation
 */
#define GG_DECODER  0
#define OA_DECODER  1
#define BD_DECODER  2
#define CBD_DECODER 3

/*
 * Type of scheduling algorithms for SNC recoding
 * TRIV - Schedule with trivally randomly (even schedule an empty buffer)
 * RAND - Randomly schedule generations with non-empty buffer
 * MLPI - Maximum Local Potential Innovativeness scheduling
 */
#define TRIV_SCHED    0
#define RAND_SCHED    1
#define MLPI_SCHED    2

struct snc_context;     // Sparse network code encode context

struct snc_packet {
    int         gid;    // subgeneration id;
    GF_ELEMENT  *coes;  // SIZE_G coding coefficients of coded packet
    GF_ELEMENT  *syms;  // SIZE_P symbols of coded packet
};

struct snc_parameter {
    long    datasize;
    double  pcrate;
    int     size_b;
    int     size_g;
    int     size_p;
    int     type;
    int     bpc;         // binary precode
    int     bnc;         // binary network coding
    int     sys;         // systematic code
};

// Metainfo of the data to be snc-coded
struct snc_metainfo {
    long    datasize;   // Data size in bytes.
    double  pcrate;     // precode rate
    int     size_b;
    int     size_g;
    int     size_p;
    int     type;       // Code type
    int     bpc;        // binary precode
    int     bnc;        // binary network coding
    int     sys;        // systematic code
    int     snum;       // Number of source packets splitted from the data.
    int     cnum;       // Number of parity-check packets (cnum ~= snum * pcrate)
    int     gnum;       // Number of subgenerations
};

struct snc_decoder;     // Sparse network code decoder

struct snc_buffer;      // Buffer for storing snc packets

/*------------------------------- sncEncoder -------------------------------*/
/**
 * Create encode context from a message buffer pointed by buf. Code parameters
 * are provided in sp.
 *
 * Return Values:
 *   On success, a pointer to the allocated encode context is returned;
 *   On error, NULL is returned, and errno is set appropriately.
 **/
struct snc_context *snc_create_enc_context(unsigned char *buf, struct snc_parameter sp);

// Get code metainfo of an encode context
struct snc_metainfo *snc_get_metainfo(struct snc_context *sc);

// Load file content into encode context
int snc_load_file_to_context(const char *filepath, long start, struct snc_context *sc);

// Free up encode context
void snc_free_enc_context(struct snc_context *sc);

// Restore data in the encode context to a char buffer
unsigned char *snc_recover_data(struct snc_context *sc);

// Free the buffer of the recovered data
void snc_free_recovered(unsigned char *data);

// Restore data in the encode context to a file (append if file exists)
long snc_recover_to_file(const char *filepath, struct snc_context *sc);

// Allocate an snc packet with coes and syms being zero
struct snc_packet *snc_alloc_empty_packet(struct snc_metainfo *meta);

// Generate an snc packet from the encode context
struct snc_packet *snc_generate_packet(struct snc_context *sc);

// Generate an snc packet to the memory of an existing snc_packet struct
int snc_generate_packet_im(struct snc_context *sc, struct snc_packet *pkt);

// Free up an snc packet
void snc_free_packet(struct snc_packet *pkt);

// Print encode/decode summary of an snc (for benchmarking)
void print_code_summary(struct snc_context *sc, int overhead, long long operations);

/*------------------------------- sncDecoder -------------------------------*/
/**
 * Create an snc decoder given code parameter and decoder type
 *   4 decoders are supported:
 *      GG_DECODER
 *      OA_DECODER
 *      BD_DECODER
 *      CBD_DECODER
 */
struct snc_decoder *snc_create_decoder(struct snc_parameter sp, int d_type);

// Get the encode context that the decoder is working on/finished.
struct snc_context *snc_get_enc_context(struct snc_decoder *decoder);

// Feed decoder with an snc packet
void snc_process_packet(struct snc_decoder *decoder, struct snc_packet *pkt);

// Check whether the decoder is finished
int snc_decoder_finished(struct snc_decoder *decoder);

// Overhead of code
// Returns the number of received packets of the decoder
int snc_code_overhead(struct snc_decoder *decoder);

// Decode cost of the decoder
// Returns the number of finite field operations the decoder has performed
long long snc_decode_cost(struct snc_decoder *decoder);

// Free decoder memory
void snc_free_decoder(struct snc_decoder *decoder);

// Save decoder context into a file
long snc_save_decoder_context(struct snc_decoder *decoder, const char *filepath);

// Restore decoder from the context stored in a file
struct snc_decoder *snc_restore_decoder(const char *filepath);

/*----------------------------- sncRecoder ------------------------------*/
/**
 * Create a buffer for storing snc packets.
 *   Arguments:
 *     snc_metainfo - meta info of the snc code
 *     bufsize      - buffer size of each subgeneration
 *   Return Value:
 *     Pointer to the buffer on success; NULL on error
 **/
struct snc_buffer *snc_create_buffer(struct snc_metainfo meta, int bufsize);

// Save an snc packet to an snc buffer
void snc_buffer_packet(struct snc_buffer *buffer, struct snc_packet *pkt);

// Recode an snc packet from an snc buffer
struct snc_packet *snc_recode_packet(struct snc_buffer *buffer, int sched_t);

// Free snc buffer
void snc_free_buffer(struct snc_buffer *buffer);

#endif /* SNC_H */
