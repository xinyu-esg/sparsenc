#include <stdlib.h>
#include <string.h>
#if defined(INTEL_SSSE3)
#include <tmmintrin.h>
#endif
#include "galois.h"
static int constructed = 0;
static uint8_t galois_log_table[1<<GF_POWER];
static uint8_t galois_ilog_table[(1<<GF_POWER)];
static uint8_t galois_mult_table[(1<<GF_POWER)*(1<<GF_POWER)];
static uint8_t galois_divi_table[(1<<GF_POWER)*(1<<GF_POWER)];

/* Two half tables are used for SSE multiply_add_region*/
#if defined(INTEL_SSSE3)
static uint8_t galois_half_mult_table_high[(1<<GF_POWER)][(1<<(GF_POWER/2))];
static uint8_t galois_half_mult_table_low[(1<<GF_POWER)][(1<<(GF_POWER/2))];
#endif

static int primitive_poly_1  = 1;
static int primitive_poly_2  = 07;
static int primitive_poly_4  = 023;
static int primitive_poly_8  = 0435;    /* 100 011 101: x^8 + x^4 + x^3 + x^2 + 1 */
static int galois_create_log_table(int npower);
static int galois_create_mult_table(int npower);

int GFConstructed() {
	return constructed;
}

int constructField(int npower) 
{
	if (constructed)
		return 0;
	else {
		if (galois_create_mult_table(npower) < 0) {
			perror("constructField");
        	exit(1);
    	}
#if defined(INTEL_SSSE3)
		/* 
		 * Create half tables for SSE multiply_add_region. Currently only support GF(256).
		 */
		if (npower == 8) {
			int a, b, c, d;
			int pp = primitive_poly_8;
  			for (a = 1; a < (1<<(npower/2)) ; a++) {
    			b = 1;
    			c = a;
    			d = (a << (npower/2));
    			do {
      				galois_half_mult_table_low[b][a] = c;
					galois_half_mult_table_high[b][a] = d;
      				b <<= 1;
      				if (b & (1<<npower)) b ^= pp;
      					c <<= 1;
      				if (c & (1<<npower)) c ^= pp;
      					d <<= 1;
      				if (d & (1<<npower)) d ^= pp;
    			} while (c != a);
  			}	
		}
#endif
		constructed = 1;
	}
	return 0;
}

static int galois_create_log_table(int npower)
{
  	int j, b;
	int nw	 = (1<<GF_POWER);
	int nwml = (1<<GF_POWER)-1;
	int m = npower;
	int gf_poly;

	if (m == 1) {
		gf_poly = primitive_poly_1;
		nw		=  1 << 1;
		nwml 	= (1 << 1) - 1;
	} else if (m == 2) {
		gf_poly = primitive_poly_2;
		nw		=  1 << 2;
		nwml 	= (1 << 2) - 1;
	} else if (m == 4) {
		gf_poly = primitive_poly_4;
		nw		=  1 << 4;
		nwml	= (1 << 4) - 1;
	} else if (m == 8) {
		gf_poly = primitive_poly_8;
		nw		=  1 << 8;
		nwml	= (1 << 8) - 1;
	} else {
		fprintf(stderr, "ERROR! We only accept field size 2^[1, 2, 4, 8, 16]\n");
		return 0;
	}


 	for (j=0; j<nw; j++) {
    	galois_log_table[j] = nwml;
    	galois_ilog_table[j] = 0;
  	} 
  
  	b = 1;
  	for (j=0; j<nwml; j++) {
    	if (galois_log_table[b] != nwml) {
      		fprintf(stderr, "Galois_create_log_tables Error: j=%d, b=%d, B->J[b]=%d, J->B[j]=%d (0%o)\n", j, b, galois_log_table[b], galois_ilog_table[j], (b << 1) ^ gf_poly);
      		exit(1);
    	}
    	galois_log_table[b] = j;
    	galois_ilog_table[j] = b;
    	b = b << 1;
    	if (b & nw) 
			b = (b ^ gf_poly) & nwml;
  	}
  
  	return 0;
}

static int galois_create_mult_table(int npower)
{
	int j, x, y, logx;
	int nw = (1<<GF_POWER);

	// create tables
    if (galois_create_log_table(npower) < 0) {
		fprintf(stderr, "create log/ilog tables failed\n");
		return -1;
	}
 
	/* Set mult/div tables for x = 0 */
  	j = 0;
  	galois_mult_table[j] = 0;   /* y = 0 */
  	galois_divi_table[j] = -1;
  	j++;
  	for (y=1; y<nw; y++) {   /* y > 0 */
    	galois_mult_table[j] = 0;
    	galois_divi_table[j] = 0;
    	j++;
  	}
  
  	for (x=1; x<nw; x++) {  /* x > 0 */
    	galois_mult_table[j] = 0; /* y = 0 */
    	galois_divi_table[j] = -1;
    	j++;
    	logx = galois_log_table[x];
    
		for (y=1; y<nw; y++) {  /* y > 0 */
			int tmp;
			tmp = logx + galois_log_table[y];
			if (tmp >= ((1<<GF_POWER) - 1))
				tmp -= ((1<<GF_POWER) - 1);				// avoid cross the boundary of log/ilog tables
 			galois_mult_table[j] = galois_ilog_table[tmp];

			tmp = logx - galois_log_table[y];
			while (tmp < 0)
				tmp += ((1<<GF_POWER) - 1);
      		galois_divi_table[j] = galois_ilog_table[tmp]; 
      		
			j++;
    	}
  	}

  	return 0;
}

// add operation over GF(2^m)
inline uint8_t galois_add(uint8_t a, uint8_t b)
{
	return a ^ b;
}

inline uint8_t galois_sub(uint8_t a, uint8_t b)
{
	return a ^ b;
}

inline uint8_t galois_multiply(uint8_t a, uint8_t b, int npower)
{
	if (a ==0 || b== 0)
		return 0;

	if (a == 1)
		return b;
	else if (b == 1)
		return a;
	
	uint8_t result = galois_mult_table[(a<<npower) | b];
	return result;
}

// return a/b
inline uint8_t galois_divide(uint8_t a, uint8_t b, int npower)
{
	if (b == 0) {
		fprintf(stderr, "ERROR! Divide by ZERO!\n");
		return -1;
	}
	
	if (a == 0)
		return 0;

	if (b == 1)
		return a;

	uint8_t result =  galois_divi_table[(a<<npower) | b];
	return result;
}

/*
 * When SSE is enabled, use SSE instructions to do multiply_add_region
 */
void galois_multiply_add_region(uint8_t *dst, uint8_t *src, uint8_t multiplier, int bytes, int npower)
{
	if (multiplier == 0) {
		// add nothing to bytes starting from *dst, just return
		return;
	}
	int i;
#if	defined(INTEL_SSSE3)
  	uint8_t *sptr, *dptr, *top;
	sptr = src;
	dptr = dst;
	top  = src + bytes;

  	uint8_t *bh, *bl;
	__m128i mth, mtl, loset;
	if (multiplier != 1) {
		/* half tables only needed for multiplier != 1 */
		bh = (uint8_t*) galois_half_mult_table_high;
		bh += (multiplier << 4);
  		bl = (uint8_t*) galois_half_mult_table_low;
  		bl += (multiplier << 4);
		// read split tables as 128-bit values
  		mth = _mm_loadu_si128((__m128i *)(bh));
  		mtl = _mm_loadu_si128((__m128i *)(bl));
  		loset = _mm_set1_epi8(0x0f);
	}

	__m128i va, vb, r, t1;
    while (sptr < top) 
	{
		if (sptr + 16 > top) {
			/* remaining data doesn't fit into __m128i, do not use SSE */
			for (i=0; i<top-sptr; i++) {
				if (multiplier == 1)
					*(dptr+i) ^= *(sptr+i);
				else
					*(dptr+i) ^= galois_mult_table[((*(sptr+i))<<npower) | multiplier];
			}
			break;
		}
		va = _mm_loadu_si128 ((__m128i *)(sptr));
		if (multiplier == 1) {
			/* just XOR */
			vb = _mm_loadu_si128 ((__m128i *)(dptr));
			vb = _mm_xor_si128(va, vb);
      		_mm_storeu_si128 ((__m128i *)(dptr), vb);
		} else {
			/* use half tables */
      		t1 = _mm_and_si128 (loset, va);
      		r = _mm_shuffle_epi8 (mtl, t1);
      		va = _mm_srli_epi64 (va, 4);
      		t1 = _mm_and_si128 (loset, va);
      		r = _mm_xor_si128 (r, _mm_shuffle_epi8 (mth, t1));
      		va = _mm_loadu_si128 ((__m128i *)(dptr));
      		r = _mm_xor_si128 (r, va);
      		_mm_storeu_si128 ((__m128i *)(dptr), r);
		}
      	dptr += 16;
      	sptr += 16;
    }
	return;
#else
	if (multiplier == 1) {
		for (i=0; i<bytes; i++)
			dst[i] ^= src[i];
		return;
	}

    for (i = 0; i < bytes; i++) 
      	dst[i] ^= galois_mult_table[(src[i]<<npower) | multiplier];
	return;
#endif
}

/*
 * Muliply a region of elements with multiplier. When SSE is available, use it
 */
void galois_multiply_region(uint8_t *src, uint8_t multiplier, int bytes, int npower)
{
	if (multiplier == 0) {
		memset(src, 0, sizeof(uint8_t)*bytes);
		return;
	} else if (multiplier == 1) {
		return;
	}
#if defined(INTEL_SSSE3)
  	uint8_t *sptr, *top;
	sptr = src;
	top  = src + bytes;

  	uint8_t *bh, *bl;
	/* half tables only needed for multiplier != 1 */
	bh = (uint8_t*) galois_half_mult_table_high;
	bh += (multiplier << 4);
  	bl = (uint8_t*) galois_half_mult_table_low;
  	bl += (multiplier << 4);
	// read split tables as 128-bit values
  	__m128i mth = _mm_loadu_si128((__m128i *)(bh));
  	__m128i mtl = _mm_loadu_si128((__m128i *)(bl));
  	__m128i loset = _mm_set1_epi8(0x0f);

	__m128i va, r, t1;
    while (sptr < top) 
	{
		if (sptr + 16 > top) {
			/* remaining data doesn't fit into __m128i, do not use SSE */
			for (int i=0; i<top-sptr; i++)
				*(sptr+i) = galois_mult_table[((*(sptr+i))<<npower) | multiplier];
			break;
		}
		va = _mm_loadu_si128 ((__m128i *)(sptr));
      	t1 = _mm_and_si128 (loset, va);
      	r = _mm_shuffle_epi8 (mtl, t1);
      	va = _mm_srli_epi64 (va, 4);
      	t1 = _mm_and_si128 (loset, va);
      	r = _mm_xor_si128 (r, _mm_shuffle_epi8 (mth, t1));
		_mm_storeu_si128 ((__m128i *)(sptr), r);
      	sptr += 16;
    }
	return;
#else
	for (int i=0; i<bytes; i++) 
		src[i] = galois_mult_table[((src[i])<<npower) | multiplier];
	return;
#endif
}
