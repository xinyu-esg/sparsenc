#include <stdio.h>
#include <string.h>
#include "gncEncoder.h"

void print_code_summary(struct gnc_metainfo *meta, int overhead, long long operations)
{
	char typestr[20];
	switch(meta->type) {
		case RAND_GNC_CODE:
			strcpy(typestr, "RAND");
			break;
		case BAND_GNC_CODE:
			strcpy(typestr, "BAND");
			break;
		default:
			strcpy(typestr, "UNKNOWN");
	}
	printf("datasize: %d ", meta->datasize);
	printf("snum: %d ", meta->snum);
	printf("size_b: %d ", meta->size_b);
	printf("size_g: %d ", meta->size_g);
	printf("size_p: %d ", meta->size_p);
	printf("type: %s ", typestr);
	printf("overhead: %.3f ", (double) overhead/meta->snum);
	printf("computation: %f\n", (double) operations/meta->snum/meta->size_p);
}
