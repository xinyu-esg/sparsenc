/*----------------------- bipartite.c ----------------------
 *
 *  Bipartite graph functions for LDPC precoding.
 *  We use LDPC code specified in Raptor Codes Standard
 *  Implementation. Refer to:
 *    1, Sec 5.4.2.3 in RFC5053 "Raptor Forward Error
 *       Correction Scheme for Object Delivery" by Luby et. al.
 *    2, Sec 3.2.3 in Foundations and Trends in Comm. and
 *       Info. Theory "Raptor Codes" by Shokrollahi et. al.
 *  for detailed information.
 *
 *----------------------------------------------------------*/
#include "common.h"
#include <math.h>
static int is_prime(int number);
static int include_left_node(int l_index, int r_index, BP_graph *graph);

// The number of required LDPC check symbols given the number of source packets.
int number_of_checks(int snum, double r)
{
    if (r == 0)
        return (0);
    int x = (int) floor( sqrt( 2 * snum ) );
    while ( x * (x - 1) < 2 * snum )
        x++;

    int c = (int) ceil( r * snum ) + x; // In Raptor code r=0.01
    while ( !is_prime(c) )
        c++;

    return c;
}

// construct LDPC graph
int create_bipartite_graph(BP_graph *graph, int nleft, int nright)
{
    int LDPC_SYS = nleft;
    int S        = nright;
    if (S == 0)
        return 0;

    int i, j;
    graph->nleft  = nleft;
    graph->nright = nright;
    if ((graph->l_nbrs_of_r = calloc(nright, sizeof(NBR_nodes*))) == NULL)
        goto failure;
    for (i=0; i<nright; i++) {
        if ((graph->l_nbrs_of_r[i] = malloc(sizeof(NBR_nodes))) == NULL)
            goto failure;
        graph->l_nbrs_of_r[i]->first = graph->l_nbrs_of_r[i]->last = NULL;
    }

    if ((graph->r_nbrs_of_l = calloc(nleft, sizeof(NBR_nodes*))) == NULL)
        goto failure;
    for (i=0; i<nleft; i++) {
        if ((graph->r_nbrs_of_l[i] = malloc(sizeof(NBR_nodes))) == NULL)
            goto failure;
        graph->r_nbrs_of_l[i]->first = graph->r_nbrs_of_l[i]->last = NULL;
    }

    char *hdpc = getenv("SNC_PRECODE");
    if (hdpc != NULL && strcmp(hdpc, "HDPC") == 0) {
        // A reference bipartite graph, which is much highly dense. This is only used 
        // when SNC_PRECODE env var is set to HDPC. The env var is ONLY for development
        // and testing use.
        for (i=0; i<S; i++) {
            for (j=0; j<LDPC_SYS; j++) {
                int included = 1;
                if (graph->binaryce == 1) {
                    if (snc_rand() % 2 == 0)
                        included = 0;
                } else {
                    if (snc_rand() % 256 == 0)
                        included = 0;
                }
                if (included) {
                    if (include_left_node(j, i, graph) < 0)
                        goto failure;
                }
            }
        }
        return 0;
    }

    // By default use circulant LDPC code which is used by Raptor code
    int a, b;
    int touching_edge = 0;
    for (i=0; i<ceil((double) LDPC_SYS/S); i++) {
        if (touching_edge == 1)
            break;

        // assign non-zero positions for the first column in each circulant matrix
        // each check node connects to exactly 3 left nodes
        // 1, P[0][i*S] = 1;
        if (include_left_node(i*S, 0, graph) < 0)
            goto failure;
        //2, P[a-1][i*S] = 1;
        a = (((i+1)+1)%S == 0) ? S : ((i+1)+1)%S;
        if (include_left_node(i*S, a-1, graph) < 0)
            goto failure;
        //3, P[b-1][i*S] = 1;
        b = ((2*(i+1)+1)%S == 0) ? S : (2*(i+1)+1)%S;
        if (include_left_node(i*S, b-1, graph) < 0)
            goto failure;

        // circulant part
        for (j=1; j<S; j++) {
            if (i*S + j >= LDPC_SYS) {
                touching_edge = 1;
                break;
            }
            // shift down the non-zero positions of previous columns in the circulant matrix
            //1, P[0+j][i*S+j] = 1;
            if (include_left_node(i*S+j, 0+j, graph) < 0)
                goto failure;
            //2, P[a-1][i*S+j] = 1;
            a = (((i+1)+1+j)%S == 0) ? S : ((i+1)+1+j)%S;
            if (include_left_node(i*S+j, a-1, graph) < 0)
                goto failure;
            //3, P[b-1][i*S+j] = 1;
            b = ((2*(i+1)+1+j)%S == 0) ? S : (2*(i+1)+1+j)%S;
            if (include_left_node(i*S+j, b-1, graph) < 0)
                goto failure;
        }
    }
    return 0;

failure:
    free_bipartite_graph(graph);
    return -1;
}

// include left node index in the LDPC graph
static int include_left_node(int l_index, int r_index, BP_graph *graph)
{
    // Coding coefficient associated with the edge
    GF_ELEMENT ce;
    if (graph->binaryce == 1) {
        ce = 1;
    } else {
        ce = (GF_ELEMENT) (snc_rand() % 255 + 1); // Value range: [1-255]
    }
    // Record neighbor of a right-side node
    NBR_node *nb = calloc(1, sizeof(NBR_node));
    if (nb == NULL)
        return -1;
    nb->data = l_index;
    nb->ce   = ce;
    nb->next = NULL;
    append_to_list(graph->l_nbrs_of_r[r_index], nb);

    // Also neighbour of the left-side node
    NBR_node *check_nb = calloc(1, sizeof(NBR_node));
    if (nb == NULL)
        return -1;
    check_nb->data = r_index;
    check_nb->ce   = ce;
    check_nb->next = NULL;
    append_to_list(graph->r_nbrs_of_l[l_index], check_nb);
    return 0;
}

void free_bipartite_graph(BP_graph *graph)
{
    if (graph == NULL)
        return;
    int i;
    if (graph->r_nbrs_of_l != NULL) {
        for (i=0; i<graph->nleft; i++)
            free_list(graph->r_nbrs_of_l[i]);
        free(graph->r_nbrs_of_l);
    }
    if (graph->l_nbrs_of_r != NULL) {
        for (i=0; i<graph->nright; i++)
            free_list(graph->l_nbrs_of_r[i]);
        free(graph->l_nbrs_of_r);
    }
    free(graph);
}

static int is_prime(int number)
{
    int i;
    for (i=2; i*i<=number; i++) {
        if (number % i == 0)
            return 0;
    }
    return 1;
}
