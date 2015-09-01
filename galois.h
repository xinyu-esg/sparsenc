#ifndef _CCFD_GALOIS_H
#define _CCFD_GALOIS_H
#define GF_ORDER	8			/* Order of Galois field, don't change it */
typedef unsigned char GF_ELEMENT;
// Galois field arithmetic routines
int constructField(int order);
GF_ELEMENT galois_add(GF_ELEMENT a, GF_ELEMENT b);
GF_ELEMENT galois_sub(GF_ELEMENT a, GF_ELEMENT b);
GF_ELEMENT galois_multiply(GF_ELEMENT a, GF_ELEMENT b, int order);
GF_ELEMENT galois_divide(GF_ELEMENT a, GF_ELEMENT b, int order);
void galois_multiply_add_region(GF_ELEMENT *dst, GF_ELEMENT *src, GF_ELEMENT multiplier, int bytes, int order);
#endif
