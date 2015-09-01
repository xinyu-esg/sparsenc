#include <stdio.h>
#include <stdlib.h>
#include "galois.h"
static int constructed = 0;
static GF_ELEMENT galois_log_table[1<<GF_ORDER];
static GF_ELEMENT galois_ilog_table[(1<<GF_ORDER)];
static GF_ELEMENT galois_mult_table[(1<<GF_ORDER)*(1<<GF_ORDER)];
static GF_ELEMENT galois_divi_table[(1<<GF_ORDER)*(1<<GF_ORDER)];

static int primitive_poly_1  = 1;
static int primitive_poly_2  = 07;
static int primitive_poly_4  = 023;
static int primitive_poly_8  = 0435;
static int galois_create_log_table(int order);
static int galois_create_mult_table(int order);

int GFConstructed() {
	return constructed;
}

int constructField(int order) //Construct the Galois field, that is, construct the two tables
{
	static int constructed = 0;			// ensure Galois field is constructed only once
	if (constructed)
		return 0;
	else {
		if (galois_create_mult_table(order) < 0) {
			perror("constructField");
        	exit(1);
    	}
		constructed = 1;
	}
	return 0;
}

static int galois_create_log_table(int order)
{
  	int j, b;
	int nw	 = (1<<GF_ORDER);
	int nwml = (1<<GF_ORDER)-1;
	int m = order;
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
		printf("ERROR! We only accept field size 2^[1, 2, 4, 8, 16]\n");
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

static int galois_create_mult_table(int order)
{
	int j, x, y, logx;
	int nw = (1<<GF_ORDER);

	// create tables
    if (galois_create_log_table(order) < 0) {
		printf("create log/ilog tables failed\n");
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
			if (tmp >= ((1<<GF_ORDER) - 1))
				tmp -= ((1<<GF_ORDER) - 1);				// avoid cross the boundary of log/ilog tables
 			galois_mult_table[j] = galois_ilog_table[tmp];

			tmp = logx - galois_log_table[y];
			while (tmp < 0)
				tmp += ((1<<GF_ORDER) - 1);
      		galois_divi_table[j] = galois_ilog_table[tmp]; 
      		
			j++;
    	}
  	}

  	return 0;
}



// add operation over GF(2^m)
GF_ELEMENT galois_add(GF_ELEMENT a, GF_ELEMENT b)
{
	return a ^ b;
}

GF_ELEMENT galois_sub(GF_ELEMENT a, GF_ELEMENT b)
{
	return a ^ b;
}

GF_ELEMENT galois_multiply(GF_ELEMENT a, GF_ELEMENT b, int order)
{
	if (a ==0 || b== 0)
		return 0;
	
	GF_ELEMENT result = galois_mult_table[(a<<order) | b];
	return result;
}

// return a/b
GF_ELEMENT galois_divide(GF_ELEMENT a, GF_ELEMENT b, int order)
{
	if (b == 0) {
		printf("ERROR! Divide by ZERO!\n");
		return -1;
	}
	
	if (a == 0)
		return 0;

	GF_ELEMENT result =  galois_divi_table[(a<<order) | b];
	return result;
}

void galois_multiply_add_region(GF_ELEMENT *dst, GF_ELEMENT *src, GF_ELEMENT multiplier, int bytes, int order)
{
	int i;
	if (multiplier == 0) {
		// add nothing to bytes starting from *dst, just return
		return;
	}
	if (multiplier == 1) {
		for (i=0; i<bytes; i++)
			dst[i] ^= src[i];
		return;
	}

    for (i = 0; i < bytes; i++) 
      dst[i] ^= galois_mult_table[(src[i]<<order) | multiplier];
}
