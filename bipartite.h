#ifndef _CCFD_BIPARTITE_H
#define _CCFD_BIPARTITE_H
#include "common.h"
typedef struct node      NBR_node;
typedef struct node_list NBR_nodes;

// Bipartitle graph for LDPC code
typedef struct bipartite_graph {
	int		  nleft;
	int		  nright;
	NBR_nodes **l_nbrs_of_r;							// left side neighbours of right
	NBR_nodes **r_nbrs_of_l;							// right side neighbours of left					
} BP_graph;

void create_bipartite_graph(BP_graph *graph, int nleft, int nright);
void free_bipartite_graph(BP_graph *graph);
#endif /* _CCFD_BIPARTITE_H */
