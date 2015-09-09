#include "common.h"
#include "bipartite.h"
#include <math.h>
static void include_left_node(int l_index, int r_index, BP_graph *graph);
// The follow functions add check packets to source packets. Check packets are constructed using LDPC
// precoding. We use LDPC precoder specified in Raptor Codes Standard Implementation
// refer to: 
// 1, Sec 5.4.2.3 in RFC5053 "Raptor Forward Error Correction Scheme for Object Delivery" by Luby et. al.
// 2, Sec 3.2.3 in Foundations and Trends in Comm. and Info. Theory "Raptor Codes" by Shokrollahi et. al. 
// construct LDPC graph, using LDPC_SYS, S macros
void create_bipartite_graph(BP_graph *graph, int nleft, int nright)
{
	int LDPC_SYS = nleft;
	int S 		 = nright;
	if (S == 0)
		return;

	int i, j;
	graph->nleft  = nleft;
	graph->nright = nright;
	graph->l_nbrs_of_r = malloc(sizeof(NBR_nodes*) * nright);
	for (i=0; i<nright; i++) {
		graph->l_nbrs_of_r[i] = malloc(sizeof(NBR_nodes));
		graph->l_nbrs_of_r[i]->first = graph->l_nbrs_of_r[i]->last = NULL;
	}
	graph->r_nbrs_of_l = malloc(sizeof(NBR_nodes*) * nleft);
	for (i=0; i<nleft; i++) {
		graph->r_nbrs_of_l[i] = malloc(sizeof(NBR_nodes));
		graph->r_nbrs_of_l[i]->first = graph->r_nbrs_of_l[i]->last = NULL;
	}

	int a, b;
	int touching_edge = 0;
	for (i=0; i<ceil((double) LDPC_SYS/S); i++) {
		if (touching_edge == 1)
			break;
		
		// assign non-zero positions for the first column in each circulant matrix
		// each check node connects to exactly 3 left nodes
		// 1, P[0][i*S] = 1;
		include_left_node(i*S, 0, graph);
		//2, P[a-1][i*S] = 1;
		a = (((i+1)+1)%S == 0) ? S : ((i+1)+1)%S;
		include_left_node(i*S, a-1, graph);
		//3, P[b-1][i*S] = 1;
		b = ((2*(i+1)+1)%S == 0) ? S : (2*(i+1)+1)%S;
		include_left_node(i*S, b-1, graph);

		// circulant part
		for (j=1; j<S; j++) {
			if (i*S + j >= LDPC_SYS) {
				touching_edge = 1;
				break;
			}
			// shift down the non-zero positions of previous columns in the circulant matrix
			//1, P[0+j][i*S+j] = 1;
			include_left_node(i*S+j, 0+j, graph);
			//2, P[a-1][i*S+j] = 1;
			a = (((i+1)+1+j)%S == 0) ? S : ((i+1)+1+j)%S;
			include_left_node(i*S+j, a-1, graph);
			//3, P[b-1][i*S+j] = 1;
			b = ((2*(i+1)+1+j)%S == 0) ? S : (2*(i+1)+1+j)%S;
			include_left_node(i*S+j, b-1, graph);
		}
	}
	
}

// include left node index in the LDPC graph
static void include_left_node(int l_index, int r_index, BP_graph *graph)
{
	// Record neighbor of a right-side node
	NBR_node *nb = malloc(sizeof(NBR_node));
	nb->data = l_index;
	nb->next = NULL;
	append_to_list(graph->l_nbrs_of_r[r_index], nb);


	// Also neighbour of the left-side node
	NBR_node *check_nb = malloc(sizeof(NBR_node));
	check_nb->data = r_index;
	check_nb->next = NULL;
	append_to_list(graph->r_nbrs_of_l[l_index], check_nb);
}

void free_bipartite_graph(BP_graph *graph)
{
	int i;
	for (i=0; i<graph->nleft; i++)
		free_list(graph->r_nbrs_of_l[i]);
	for (i=0; i<graph->nright; i++) 
		free_list(graph->l_nbrs_of_r[i]);
	free(graph);
}


