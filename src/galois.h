/*-------------------- galois.h ------------------------------
 *
 * Internal header file of Galois field implementation
 *
 *-----------------------------------------------------------*/
#ifndef GALOIS_H
#define GALOIS_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define GF_POWER	8			/* Order of Galois field: 2^8, don't change it */
#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
// Galois field arithmetic routines
int constructField(int npower);
uint8_t galois_add(uint8_t a, uint8_t b);
uint8_t galois_sub(uint8_t a, uint8_t b);
uint8_t galois_multiply(uint8_t a, uint8_t b, int npower);
uint8_t galois_divide(uint8_t a, uint8_t b, int npower);
void galois_multiply_region(uint8_t *src, uint8_t multiplier, int bytes, int npower);
void galois_multiply_add_region(uint8_t *dst, uint8_t *src, uint8_t multiplier, int bytes, int npower);
#endif
