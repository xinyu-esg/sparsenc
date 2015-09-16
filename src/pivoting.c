/********************************************************************
 *
 * pivoting.c
 *
 * This library provides routines to pivot matrices of a linear 
 * system of equations. Solving the system would require much lower
 * computational cost after pivoting. Two pivoting algorithms, namely
 * inactivation and Zlatev pivoting, are implemented. One interface 
 * is exposed to caller.
 *
 ********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "galois.h"
#define ZLATEVS			3    // Number of searching rows using Zlatev's strategy 
#define IA_INIT			0    // Initial number of inactivated columns
#define IA_STEP			1    // Gradually inactivate more columns
/*
 * Pivoting algorithms use double-linked lists to store numbers of 
 * nonzeros of rows and columns of a matrix.
 */
typedef struct subscript				Subscript;
typedef struct subscripts				ssList;
struct subscript
{
	int index;
	int nonzeros;
	struct subscript *next;
	struct subscript *prev;
};

struct subscripts
{
	struct subscript *ssFirst;
	struct subscript *ssLast;
};



/*
 * Procedures to pivot matrix.
 */
static int inactivation_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots);
static int zlatev_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots);

/*
 * Reshape matrix A and B of Ax=B according to pivot sequence in (RowPivots, ColPivots)
 */
static void reshape_matrix(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, ssList *RowPivots, ssList *ColPivots);

/*
 * Reorder columns of matrix according to sequence in ColPivots
 */
static void permute_matrix_columns(int nrow, int ncolA, GF_ELEMENT **A, ssList *ColPivots);

/* 
 * Helper functions for pivoting
 * We use double-linked lists when pivoting a matrix.
 */
static void insertSubAtBeginning(ssList *sub_list, Subscript *sub);
static void insertSubAtEnd(ssList *sub_list, Subscript *sub);
static void removeSubscript(ssList *sub_list, Subscript *sub);
static void free_subscriptList(ssList *sub_list);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);

/**********************************************************************************
 * pivot_matrix_x()
 * Main interfaces to the library.
 * Input:
 * 	Ax = B
 * 	It is required that 2-D matrices A and B are stored in form of double pointers
 * 	
 * 	A ->A[0] A[0][1] A[0][2] ...
 * 		A[1] ...
 * 		.
 * 		.
 * Parameters:
 * 	nrow  - number of rows of A and B
 * 	ncolA - number of columns of A
 * 	ncolB - number of columns of B
 *
 * Return: 
 *  number of Galois field operations consumed
 *
 * Return as arguments:
 *  otoc        - An array containing the current column index of the original index
 *  	          i.e., otoc[i] = "current column index of the original i-th column"
 *  inactives   - number of inactivated columns
 *
 *  [A] will be transformed into the form of:
 * 		-                 -
 * 		| x 0 0 0 0 x x x |
 * 		| 0 x 0 0 0 x x x |
 * 		| 0 0 x 0 0 x x x |
 * 		| 0 0 0 x 0 x x x |
 * 		| 0 0 0 0 x x x x |
 * 		| 0 0 0 0 0 x x x |
 * 		| 0 0 0 0 0 x x x |
 * 		| 0 0 0 0 0 x x x |
 * 		-                -
 *  and B will be processed accordingly.
 *
 **********************************************************************************/

/*
 * Pivot matrix use inactivation pivoting
 */
long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int *otoc, int *inactives)
{
	int i, j, k;
	long operations = 0;
	// First pivoting: inactivation
	ssList *row_pivots = malloc(sizeof(ssList));
	ssList *col_pivots = malloc(sizeof(ssList));
	row_pivots->ssFirst = row_pivots->ssLast = NULL;
	col_pivots->ssFirst = col_pivots->ssLast = NULL;
	
	// Make a local copy of A and B
	GF_ELEMENT **ces_matrix = calloc(nrow, sizeof(GF_ELEMENT*));
	GF_ELEMENT **msg_matrix = calloc(nrow, sizeof(GF_ELEMENT*));
	for (i=0; i<nrow; i++) {
		ces_matrix[i] = calloc(ncolA, sizeof(GF_ELEMENT));
		memcpy(ces_matrix[i], A[i], ncolA * sizeof(GF_ELEMENT));
		msg_matrix[i] = calloc(ncolB, sizeof(GF_ELEMENT));
		memcpy(msg_matrix[i], B[i], ncolB * sizeof(GF_ELEMENT));
	}

#if defined(GNCTRACE)
	double pivoting_time=0.0;
	time_t start_pivoting, stop_pivoting;
	time(&start_pivoting);
	printf("Inactivation pivoting strategy is used. IA_INIT: %d, IA_STEP: %d.\n", IA_INIT, IA_STEP);
#endif
	int ias = inactivation_pivoting(nrow, ncolA, ces_matrix, row_pivots, col_pivots);
	*inactives = ias;
#if defined(GNCTRACE)
	printf("A total of %d/%d columns are inactivated.\n", ias, ncolA);
	time(&stop_pivoting);
	pivoting_time += difftime(stop_pivoting, start_pivoting);
	printf("Time consumed in pivoting T: %.0f seconds.\n", pivoting_time);
#endif

	// Save orders of column indices after re-ordering
	// Current-to-original mapping
	int *ctoo = (int *) malloc(sizeof(int)*ncolA);
	Subscript *ss_pt = col_pivots->ssFirst;
	for (i=0; i<ncolA; i++) {
		otoc[ss_pt->index] = i;	 // Original-to-current mapping	
		ctoo[i] = ss_pt->index;	
		ss_pt = ss_pt->next;
	}

#if defined(GNCTRACE)
	printf("Start matrix re-ordering after the inactivation...\n");
	double reordering_time = 0.0;
	time_t start_reorder, stop_reorder;
	time(&start_reorder);
#endif
	/* reorder rows and columns of ces_matrix and msg_matrix after pivoting*/
	reshape_matrix(nrow, ncolA, ncolB, ces_matrix, msg_matrix, row_pivots, col_pivots);
	free_subscriptList(row_pivots);
	free_subscriptList(col_pivots);
#if defined(GNCTRACE)
	time(&stop_reorder);
	reordering_time = difftime(stop_reorder, start_reorder);
	printf("Time consumed in matrix re-ordering after pivoting T: %.0f seconds.\n", reordering_time);
	pivoting_time += reordering_time;
	//diagonalize active part
	printf("Convert the left half of T (lower triangular) to diagonal...\n");	
#endif
	long long ops1=0;
	GF_ELEMENT quotient;
	for (i=0; i<ncolA-ias; i++) {
		for (j=i+1; j<nrow; j++) {
#if defined(GNCTRACE)
			if (ces_matrix[i][i] == 0)
				printf("The diagonal element after re-ordering is nonzero.\n");
#endif
			// process the item on (j, i)
			if (ces_matrix[j][i] != 0) {
				quotient = galois_divide(ces_matrix[j][i], ces_matrix[i][i], GF_POWER);
				ops1 += 1;
				// XOR the corresponding part in the inactive part
				galois_multiply_add_region(&(ces_matrix[j][ncolA-ias]), &(ces_matrix[i][ncolA-ias]), quotient, ias, GF_POWER);
				ops1 += ias;	// This part of matrix in processing is lower triangular in part, so operations only needed in the back
				// simultaneously do the same thing on right matrix B
				galois_multiply_add_region(msg_matrix[j], msg_matrix[i], quotient, ncolB, GF_POWER);
				ops1 += ncolB;				
				ces_matrix[j][i] = 0;			// eliminate the item
			}
		}
	}
	operations += ops1;
	
	/* Perform forward substitution on the ias x ias dense inactivated matrix. */ 
	GF_ELEMENT **ces_submatrix = calloc(ias+nrow-ncolA, sizeof(GF_ELEMENT*));
	GF_ELEMENT **msg_submatrix = calloc(ias+nrow-ncolA, sizeof(GF_ELEMENT*));
	for (i=0; i<ias+nrow-ncolA; i++){
		ces_submatrix[i] = calloc(ias, sizeof(GF_ELEMENT));
		memcpy(ces_submatrix[i], &(ces_matrix[ncolA-ias+i][ncolA-ias]), ias*sizeof(GF_ELEMENT));
		msg_submatrix[i] = calloc(ncolB, sizeof(GF_ELEMENT));
		memcpy(msg_submatrix[i], msg_matrix[ncolA-ias+i], ncolB*sizeof(GF_ELEMENT));
	}

#if defined(GNCTRACE)
	printf("Perform forward substitution on the bottom-right part of GDM (size: %d x %d)...\n", ias, ias);
#endif
	long long ops = forward_substitute(ias+nrow-ncolA, ias, ncolB, ces_submatrix, msg_submatrix);
	operations += ops;

	// save the processed matrices back to A and B
	for (i=0; i<ncolA-ias; i++) {
		memcpy(A[i], ces_matrix[i], ncolA*sizeof(GF_ELEMENT));
		memcpy(B[i], msg_matrix[i], ncolB*sizeof(GF_ELEMENT));
	}

	for (i=ncolA-ias; i<ncolA; i++) {
		memcpy(&(A[i][ncolA-ias]), ces_submatrix[i-(ncolA-ias)], ias*sizeof(GF_ELEMENT));
		memcpy(B[i], msg_submatrix[i-(ncolA-ias)], ncolB*sizeof(GF_ELEMENT));
	}

	// free ces_matrix and msg_matrix
	for (i=0; i<nrow; i++) {
		free(ces_matrix[i]);
		free(msg_matrix[i]);
	}
	free(ces_matrix);
	free(msg_matrix);

	//free ces_submatrix and msg_submatrix
	for (i=0; i<ias+nrow-ncolA; i++){
		free(ces_submatrix[i]);
		free(msg_submatrix[i]);
	}
	free(ces_submatrix);
	free(msg_submatrix);

	return operations;
}

/********************************************************************************
 * 
 * Pivot matrix for two rounds. First round uses inactivation pivoting against 
 * the whole matrix. The second round performs on the relatively denser 
 * bottom-right corner corresponding to the inactivated columns. Zlatev pivoting
 * is used for the second round.
 *
 ********************************************************************************/
long pivot_matrix_tworound(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int *otoc, int *inactives)
{
	int i, j, k;
	long operations = 0;
	// First pivoting: inactivation
	ssList *row_pivots = malloc(sizeof(ssList));
	ssList *col_pivots = malloc(sizeof(ssList));
	row_pivots->ssFirst = row_pivots->ssLast = NULL;
	col_pivots->ssFirst = col_pivots->ssLast = NULL;
#if defined(GNCTRACE)
	double pivoting_time=0.0;
	time_t start_pivoting, stop_pivoting;
	printf("Start pivoting...\n");
	time(&start_pivoting);
#endif
	
	// Make a local copy of A and B
	GF_ELEMENT **ces_matrix = calloc(nrow, sizeof(GF_ELEMENT*));
	GF_ELEMENT **msg_matrix = calloc(nrow, sizeof(GF_ELEMENT*));
	for (i=0; i<nrow; i++) {
		ces_matrix[i] = calloc(ncolA, sizeof(GF_ELEMENT));
		memcpy(ces_matrix[i], A[i], ncolA * sizeof(GF_ELEMENT));
		msg_matrix[i] = calloc(ncolB, sizeof(GF_ELEMENT));
		memcpy(msg_matrix[i], B[i], ncolB * sizeof(GF_ELEMENT));
	}

#if defined(GNCTRACE)
	printf("Inactivation pivoting strategy is used. IA_INIT: %d, IA_STEP: %d.\n", IA_INIT, IA_STEP);
#endif
	int ias = inactivation_pivoting(nrow, ncolA, ces_matrix, row_pivots, col_pivots);
	*inactives = ias;
#if defined(GNCTRACE)
	printf("A total of %d/%d columns are inactivated.\n", ias, ncolA);
	time(&stop_pivoting);
	pivoting_time += difftime(stop_pivoting, start_pivoting);
	printf("Time consumed in pivoting T: %.0f seconds.\n", pivoting_time);
#endif

	// 保存re-order过后列的顺序
	// Inactivation后第一次re-ordering，记下indices的mapping
	int *ctoo = (int *) malloc(sizeof(int)*ncolA);
	Subscript *ss_pt = col_pivots->ssFirst;
	for (i=0; i<ncolA; i++) {
		otoc[ss_pt->index] = i;		// 第i个pivot对应的原来第ss_pt->index列
		ctoo[i] = ss_pt->index;	
		ss_pt = ss_pt->next;
	}

#if defined(GNCTRACE)
	printf("Start matrix re-ordering after the inactivation...\n");
	double reordering_time = 0.0;
	time_t start_reorder, stop_reorder;
	time(&start_reorder);
#endif
	/* reorder rows and columns of ces_matrix and msg_matrix after pivoting*/
	reshape_matrix(nrow, ncolA, ncolB, ces_matrix, msg_matrix, row_pivots, col_pivots);
	free_subscriptList(row_pivots);
	free_subscriptList(col_pivots);
#if defined(GNCTRACE)
	time(&stop_reorder);
	reordering_time = difftime(stop_reorder, start_reorder);
	printf("Time consumed in matrix re-ordering after pivoting T: %.0f seconds.\n", reordering_time);
	pivoting_time += reordering_time;

	int missing_pivots = 0;
	for (i=0; i<nrow; i++) {
		if (ces_matrix[i][i] == 0) {
			missing_pivots += 1;
		}
	}
	printf("There are %d pivots missing\n", missing_pivots);
	/* Check the density of the lower half of the active part */
	long long nonzeros = 0;
	for (i=ncolA-ias; i<ncolA; i++) {
		for (j=0; j<ncolA-ias; j++)
			if (ces_matrix[i][j] != 0)
				nonzeros++;
	}
	printf("Fraction of nonzeros in the lower half of the active part: %f\n", (double) nonzeros/(ncolA-ias)/ias);
	//diagonalize active part
	printf("Convert the left half of T (lower triangular) to diagonal...\n");	
#endif
	long long ops1=0;
	GF_ELEMENT quotient;
	for (i=0; i<ncolA-ias; i++) {
		for (j=i+1; j<nrow; j++) {
#if defined(GNCTRACE)
			if (ces_matrix[i][i] == 0)
				printf("The diagonal element after re-ordering is nonzero.\n");
#endif
			// process the item on (j, i)
			if (ces_matrix[j][i] != 0) {
				quotient = galois_divide(ces_matrix[j][i], ces_matrix[i][i], GF_POWER);
				ops1 += 1;
				// XOR the corresponding part in the inactive part
				galois_multiply_add_region(&(ces_matrix[j][ncolA-ias]), &(ces_matrix[i][ncolA-ias]), quotient, ias, GF_POWER);
				ops1 += ias;	// This part of matrix in processing is lower triangular in part, so operations only needed in the back
				// simultaneously do the same thing on right matrix B
				galois_multiply_add_region(msg_matrix[j], msg_matrix[i], quotient, ncolB, GF_POWER);
				ops1 += ncolB;				
				ces_matrix[j][i] = 0;			// eliminate the item
			}
		}
	}
	operations += ops1;
	
	/* Second time pivoting on the iasxias dense inactivated matrix. */ 
#if defined(GNCTRACE)
	printf("Perform another time of pivoting after partly diagonalizing T\n");
#endif
	GF_ELEMENT **ces_submatrix = calloc(ias+nrow-ncolA, sizeof(GF_ELEMENT*));
	GF_ELEMENT **msg_submatrix = calloc(ias+nrow-ncolA, sizeof(GF_ELEMENT*));
	for (i=0; i<ias+nrow-ncolA; i++){
		ces_submatrix[i] = calloc(ias, sizeof(GF_ELEMENT));
		memcpy(ces_submatrix[i], &(ces_matrix[ncolA-ias+i][ncolA-ias]), ias*sizeof(GF_ELEMENT));
		msg_submatrix[i] = calloc(ncolB, sizeof(GF_ELEMENT));
		memcpy(msg_submatrix[i], msg_matrix[ncolA-ias+i], ncolB*sizeof(GF_ELEMENT));
	}

	ssList *row_pivots_2nd = malloc(sizeof(ssList));
	ssList *col_pivots_2nd = malloc(sizeof(ssList));
	row_pivots_2nd->ssFirst = row_pivots_2nd->ssLast = NULL;
	col_pivots_2nd->ssFirst = col_pivots_2nd->ssLast = NULL;
#if defined(GNCTRACE)
	printf("Start the second-time pivoting...\n");
#endif
	int ias_2nd = zlatev_pivoting(ias+nrow-ncolA, ias, ces_submatrix, row_pivots_2nd, col_pivots_2nd);
#if defined(GNCTRACE)
	printf("A total of %d/%d columns are further inactivated.\n", ias_2nd, ias);
#endif
	
	// Update otoc_mapping & ctoo_mapping after second-time pivoting
	// Careful!! This step is critical.
	int *partial_mappings = (int *) calloc(ias, sizeof(int));

	Subscript *ss_pt_2nd = col_pivots_2nd->ssFirst;
	for (i=0; i<ias; i++) {
		partial_mappings[i] = ctoo[ncolA-ias+ss_pt_2nd->index];
		ss_pt_2nd = ss_pt_2nd->next;
	}
	memcpy(&(ctoo[ncolA-ias]), partial_mappings, sizeof(int)*ias);
	free(partial_mappings);
	for (i=0; i<ncolA; i++) {
		otoc[ctoo[i]] = i;
	}
	free(ctoo);

	// re-order the ((ias+OHS) x ias) matrix
	reshape_matrix(ias+nrow-ncolA, ias, ncolB, ces_submatrix, msg_submatrix, row_pivots_2nd, col_pivots_2nd);
	
	// re-order the columns above it accordingly
	GF_ELEMENT **BI_matrix = calloc(nrow-ias, sizeof(GF_ELEMENT*));
	for (i=0; i<nrow-ias; i++) {
		BI_matrix[i] = calloc(ias, sizeof(GF_ELEMENT));
		memcpy(BI_matrix[i], &(ces_matrix[i][ncolA-ias]), ias*sizeof(GF_ELEMENT));
	}

	permute_matrix_columns(nrow-ias, ias, BI_matrix, col_pivots_2nd);
	free_subscriptList(row_pivots_2nd);
	free_subscriptList(col_pivots_2nd);
	// copy re-ordered submatrix B_I back
	for (i=0; i<nrow-ias; i++) 
		memcpy(&(ces_matrix[i][ncolA-ias]), BI_matrix[i], ias*sizeof(GF_ELEMENT));
	// second pivoting is done

	// 3, perform forward substitution on the re-ordered matrix
#if defined(GNCTRACE)
	printf("Perform forward substitution on the bottom-right part of GDM (size: %d x %d)...\n", ias, ias);
#endif
	long long ops = forward_substitute(ias+nrow-ncolA, ias, ncolB, ces_submatrix, msg_submatrix);
	operations += ops;

	// save the processed matrices back to A and B
	for (i=0; i<ncolA-ias; i++) {
		memcpy(A[i], ces_matrix[i], ncolA*sizeof(GF_ELEMENT));
		memcpy(B[i], msg_matrix[i], ncolB*sizeof(GF_ELEMENT));
	}

	for (i=ncolA-ias; i<ncolA; i++) {
		memcpy(&(A[i][ncolA-ias]), ces_submatrix[i-(ncolA-ias)], ias*sizeof(GF_ELEMENT));
		memcpy(B[i], msg_submatrix[i-(ncolA-ias)], ncolB*sizeof(GF_ELEMENT));
	}

	// free ces_matrix and msg_matrix
	for (i=0; i<nrow; i++) {
		free(ces_matrix[i]);
		free(msg_matrix[i]);
	}
	free(ces_matrix);
	free(msg_matrix);

	//free ces_submatrix and msg_submatrix
	for (i=0; i<ias+nrow-ncolA; i++){
		free(ces_submatrix[i]);
		free(msg_submatrix[i]);
	}
	free(ces_submatrix);
	free(msg_submatrix);

	// free BI_matrix
	for (i=0; i<nrow-ias; i++) {
		free(BI_matrix[i]);
	}
	free(BI_matrix);

	return operations;
}


/*****************************************************************************
*     Zlatev pivoting
*
* A special kind of Markowitz pivoting in which pivots are selected from 3 
* candidates who have the smallest nonzeros.
******************************************************************************/
static int zlatev_pivoting(int nrow, int ncolA, GF_ELEMENT *A[], ssList *RowPivots, ssList *ColPivots)
{	
	// Nonzeros of each rows and cols
	int *row_counts = (int *) calloc(nrow, sizeof(int));
	int *col_counts = (int *) calloc(ncolA, sizeof(int));

	int i, j, k;
	int max_row1s = 0;			// Max number of nonzeros in a row
	int max_col1s = 0;			// Max number of nonzeros in a column
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncolA; j++) {
			if (A[i][j] != 0) {
				row_counts[i] += 1;
				if (row_counts[i] > max_row1s)
					max_row1s = row_counts[i];

				col_counts[j] += 1;
				if (col_counts[j] > max_col1s)
					max_col1s = col_counts[j];
			}
		}
	}
	/*************************************************************************
	 * A list of double lists. Elements of each double list contain rows/cols 
	 * indices. Rows/cols of same double-linked list have the same number of 
	 * nonzero elements. 
	 **************************************************************************/
	ssList **RowID_lists = (ssList **) malloc(sizeof(ssList*) * (max_row1s+1)); 
	ssList **ColID_lists = (ssList **) malloc(sizeof(ssList*) * (max_col1s+1));
	for (i=0; i<max_row1s+1; i++) {
		RowID_lists[i] = (ssList *) malloc(sizeof(ssList));
		RowID_lists[i]->ssFirst = RowID_lists[i]->ssLast = NULL;
	}
	for (i=0; i<max_col1s+1; i++) {
		ColID_lists[i] = (ssList *) malloc(sizeof(ssList));
		ColID_lists[i]->ssFirst = ColID_lists[i]->ssLast = NULL;
	}
	// Travel through the matrix and populate double-linked lists
	int allzero_rows = 0;
	for (i=0; i<nrow; i++) {
		int row_count = row_counts[i];
		Subscript *rowID = malloc(sizeof(Subscript));
		rowID->index = i;
		rowID->nonzeros = row_count;
		rowID->prev = rowID->next = NULL;
		insertSubAtBeginning(RowID_lists[row_count], rowID);
		if (row_count == 0)
			allzero_rows += 1;
	}
	int allzero_cols = 0;
	for (i=0; i<ncolA; i++) {
		int col_count = col_counts[i];
		Subscript *colID = malloc(sizeof(Subscript));
		colID->index = i;
		colID->nonzeros = col_count;
		colID->prev = colID->next = NULL;
		insertSubAtBeginning(ColID_lists[col_count], colID);
		if (col_count == 0)
			allzero_cols += 1;
	}

	int nonallzero_rows = nrow - allzero_rows;
	int nonallzero_cols = ncolA - allzero_cols;
#if defined(GNCTRACE)
	printf("there are %d all-zero rows and %d all-zero cols in the matrix.\n", allzero_rows, allzero_cols);
#endif

	// Markowitz pivoting using row_counts, col_counts, RowID_lists, ColID_lists
	int pivots_found = 0;
	int removed_cols = 0;
	int toallzero_cols = 0;
	int removed_rows = 0;
	int toallzero_rows = 0;
	int singletons = 0;
	while (pivots_found != ncolA) {
		//printf("%d rows removed, %d reduced to all-zero.\n", removed_rows, toallzero_rows);
		//printf("%d cols removed, %d reduced to all-zero.\n", removed_cols, toallzero_cols);

		// 计算Markowitz count
		int potential_r  = -1;
		int potential_c  = -1;
		int potential_mc = -1;				// potential Markowitz count
		int current_mc;
		int mc_minipos;				// the minimum possible Markowitz count in each row test
		int direct_found = 0;
		// search for one pivot
		Subscript *rows_inchecking, *cols_inchecking;
		int row_id, col_id;
		
		// Search rows
		int searched_rows = 0;
		for (i=1; i<=max_row1s; i++) {
			// 先按行的非零元素多少搜索
			mc_minipos = (i - 1) * (i - 1);
			if (RowID_lists[i]->ssFirst == NULL)
				continue;

			rows_inchecking = RowID_lists[i]->ssFirst;
			while (rows_inchecking != NULL && (searched_rows < ZLATEVS)) {
				searched_rows += 1;
				row_id = rows_inchecking->index;
				// 再按列的非零元素多少做test
				for (j=1; j<=max_col1s; j++) {
					cols_inchecking = ColID_lists[j]->ssFirst;
					while (cols_inchecking != NULL) {
						col_id = cols_inchecking->index;
						if (A[row_id][col_id] != 0) {
							// we have found an entry in (row_id)-th row
							current_mc = (i-1) * (j-1);
							if (current_mc == 0) {
								
								potential_r = row_id;
								potential_c = col_id;
								if (i==1)
									singletons += 1;
								//	printf("a singleton row is found.\n");
								goto found;
							} else if (potential_mc == -1 || current_mc < potential_mc) {
								potential_r = row_id;
								potential_c = col_id;
								potential_mc = current_mc;
							}
						}
						cols_inchecking = cols_inchecking->next;
					}
				}
				rows_inchecking = rows_inchecking->next;
			}
		}

		if (potential_r == -1 || potential_c == -1) {
			//printf("error: no pivot is found, the matrix is not full-rank.\n");
#if defined(GNCTRACE)
			printf("%d rows reduced to all-zero.\n", toallzero_rows);
			printf("%d cols reduced to all-zero.\n", toallzero_cols);
			printf("(partial success) %d/%d pivots were found out of %d rows.\n", pivots_found, ncolA, nrow);
#endif
			// There are row/col being reduced to all-zero, take them
			// as pivot anyway because we have no other choices
			Subscript *sub_ptt;
			sub_ptt = ColID_lists[0]->ssFirst;
			int zerocols = 0;
			while(sub_ptt != NULL) {
				zerocols += 1;
				Subscript *newCpivot0 = malloc(sizeof(Subscript));
				newCpivot0->index = sub_ptt->index;
				newCpivot0->nonzeros = sub_ptt->nonzeros;
				newCpivot0->prev = newCpivot0->next = NULL;
				insertSubAtEnd(ColPivots, newCpivot0);
				sub_ptt = sub_ptt->next;
			}

			sub_ptt = RowID_lists[0]->ssFirst;
			//for (i=0; i<(toallzero_rows+allzero_cols); i++)
			for (i=0; i<zerocols; i++) {
				Subscript *newRpivot0 = malloc(sizeof(Subscript));
				newRpivot0->index = sub_ptt->index;
				newRpivot0->nonzeros = sub_ptt->nonzeros;
				newRpivot0->prev = newRpivot0->next = NULL;
				insertSubAtEnd(RowPivots, newRpivot0);
				sub_ptt = sub_ptt->next;
			}
#if defined(GNCTRACE)
			printf("There are %d/%d singleton rows were found as pivots.\n", singletons, pivots_found);
#endif
			//return (pivots_found+toallzero_cols);
			return ncolA;
		}

found:	
		// Found a pivot, save it
		pivots_found += 1;
		int p_r = potential_r;
		int p_c = potential_c;
		// Save coordiate (p_r, p_c)
		Subscript *newRpivot = malloc(sizeof(Subscript));
		newRpivot->index = p_r;
		newRpivot->nonzeros = row_counts[p_r];
		newRpivot->prev = newRpivot->next = NULL;
		insertSubAtEnd(RowPivots, newRpivot);
		Subscript *newCpivot = malloc(sizeof(Subscript));
		newCpivot->index = p_c;
		newCpivot->nonzeros = col_counts[p_c];
		newCpivot->prev = newCpivot->next = NULL;
		insertSubAtEnd(ColPivots, newCpivot);

		// Update row_counts[], col_counts[], RowID_lists, ColID_lists
		// 1, check nonzero elements of the row: potential_r, update numbers of nonzero elements of the correspondings columns(col_counts[] and ColID_lists).
		int nzs;
		Subscript *ss_pt;
		Subscript *ss_pt_next;
		// Note: some row/col have been eliminated, so we need to traverse Subscript when updating
		// 1, Update columns
		for (i=1; i<=max_col1s; i++) {
			ss_pt = ColID_lists[i]->ssFirst;
			while (ss_pt != NULL) {
				if (A[p_r][ss_pt->index] != 0) {
					// 该subscript对象将要被更新处理，因此要记下当前遍历的位置
					ss_pt_next = ss_pt->next;
					if (ss_pt->index == p_c) {
						removeSubscript(ColID_lists[i], ss_pt);
						removed_cols += 1;
						col_counts[ss_pt->index] -= 1;
						free(ss_pt);
					} else {
						removeSubscript(ColID_lists[i], ss_pt);
						ss_pt->nonzeros -= 1;
						insertSubAtBeginning(ColID_lists[ss_pt->nonzeros], ss_pt);
						if (ss_pt->nonzeros == 0)
							toallzero_cols += 1;
						col_counts[ss_pt->index] -= 1;
					}
					ss_pt = ss_pt_next;
				} else {
					ss_pt = ss_pt->next;
				}
			}
		}
		// 2, Update rows
		for (j=1; j<=max_row1s; j++) {
			ss_pt = RowID_lists[j]->ssFirst;
			while (ss_pt != NULL) {
				if (A[ss_pt->index][p_c] != 0) {
					// 该subscript对象将要被更新处理，因此要记下当前遍历的位置
					ss_pt_next = ss_pt->next;
					if (ss_pt->index == p_r) {
						removeSubscript(RowID_lists[j], ss_pt);
						removed_rows += 1;
						row_counts[ss_pt->index] -= 1;
						free(ss_pt);
					} else {
						removeSubscript(RowID_lists[j], ss_pt);
						ss_pt->nonzeros -= 1;
						if (ss_pt->nonzeros == 0)
							toallzero_rows += 1;
						insertSubAtBeginning(RowID_lists[ss_pt->nonzeros], ss_pt);
						row_counts[ss_pt->index] -= 1;
					}
					ss_pt = ss_pt_next;
				} else {
					ss_pt = ss_pt->next;
				}
			}
		}

	}

	free(row_counts);
	free(col_counts);
	for (i=0; i<max_row1s+1; i++)
		free_subscriptList(RowID_lists[i]);
	for (i=0; i<max_col1s+1; i++)
		free_subscriptList(ColID_lists[i]);
#if defined(GNCTRACE)
	printf("%d rows reduced to all-zero.\n", toallzero_rows);
	printf("%d cols reduced to all-zero.\n", toallzero_cols);
	printf("(full success) %d/%d pivots were found out of %d rows.\n", pivots_found, ncolA, nrow);
	printf("There are %d/%d singleton rows were found as pivots.\n", singletons, pivots_found);
#endif
	return pivots_found;
}

/*********************************************************************************************************
 *      Inactivation pivoting
 * Use progressive inactivation to do pivoting. Only select pivots from singleton rows of the "residual"
 * matrix during the pivoting. If no singleton rows can be found, inactivate some columns until singletons
 * can be found:
 *  1) inactivate some columns and only perform pivoting on the rest of the "active" sub-matrix
 *  2) if singleton row cannot be found in the middle of pivoting, declare more inactive columns
 *  3) given the structure (heavier columns are in the back), declare inactive columns from the back
 *********************************************************************************************************/
static int inactivation_pivoting(int nrow, int ncolA, GF_ELEMENT *A[], ssList *RowPivots, ssList *ColPivots)
{
#if defined(GNCTRACE)
	printf("Pivoting matrix of size %d x %d via inactivation.\n", nrow, ncolA);
#endif

	int i, j, k;
	// 对矩阵中初始非零元素进行计数
	int *row_counts;
	int *col_counts;
	row_counts = (int *) malloc(nrow  * sizeof(int));
	col_counts = (int *) malloc(ncolA * sizeof(int));
	memset(row_counts, 0, nrow  * sizeof(int));
	memset(col_counts, 0, ncolA * sizeof(int));

	int max_row1s = 0;					// 记录初始矩阵里行中非零元素数目的最大值
	int max_col1s = 0;					// 记录初始矩阵里列中非零元素数目的最大值
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncolA; j++) {
			if (A[i][j] != 0) {
				row_counts[i] += 1;
				if (row_counts[i] > max_row1s)
					max_row1s = row_counts[i];

				col_counts[j] += 1;
				if (col_counts[j] > max_col1s)
					max_col1s = col_counts[j];
			}
		}
	}
	// 创建指针数组，每个指针指向一个双向链表，同一个链表中的列具有相同数目的非零元素
	// 该链表用来保存active的列的标号，按非零元素个数排列是为了方便inactivate含非零元素最多的列
	ssList **ColID_lists = (ssList **) malloc(sizeof(ssList*) * (max_col1s+1)); 
	for (i=0; i<max_col1s+1; i++) {
		ColID_lists[i] = (ssList *) malloc(sizeof(ssList));
		ColID_lists[i]->ssFirst = ColID_lists[i]->ssLast = NULL;
	}
	for (i=0; i<ncolA; i++) {
		int col_count = col_counts[i];
		Subscript *colID = malloc(sizeof(Subscript));
		colID->index = i;
		colID->nonzeros = col_count;
		colID->prev = colID->next = NULL;
		insertSubAtBeginning(ColID_lists[col_count], colID);
	}

	// three possible "states" of each column
	// 0 - active
	// 1 - inactivated
	// 2 - removed because an entry of it was chosen as the pivot
	GF_ELEMENT *inactives;				// an array record the state of columns
	inactives = (GF_ELEMENT *) malloc(sizeof(GF_ELEMENT) * ncolA);
	memset(inactives, 0, sizeof(GF_ELEMENT)*ncolA);

	int inactivated = 0;		// number of inactivated columns
	int active = ncolA;			// number of active columns
	
	int p_r;					// used to store row index of the chosen pivot 
	int p_c;					// used to store col index of the chosen pivot
	int singleton_r_found;
	int ia;
	int selected_pivots = 0;
	while (active != 0) {
		p_r = -1;
		p_c = -1;
		singleton_r_found = 0;
		for (i=0; i<nrow; i++) {
			if (row_counts[i] == 1) {
				singleton_r_found = 1;
				//printf("singleton row is found at row %d.\n", i);
				p_r = i;
				break;
			}
		}

		if (singleton_r_found == 1) {
			// a singleton row is found, store the pivot
			for (j=0; j<ncolA; j++) {
				if ((inactives[j]==0) && (A[p_r][j]!=0)) {
					p_c = j;
					break;
				}
			}
			if (p_c == -1) {
				printf("error: failed to find the nonzero element in the singlton row.\n");
				exit(1);
			}
			// save the pivot
			// 保存该pivot的坐标
			Subscript *newRpivot = malloc(sizeof(Subscript));
			newRpivot->index = p_r;
			newRpivot->nonzeros = 1;
			newRpivot->prev = newRpivot->next = NULL;
			insertSubAtEnd(RowPivots, newRpivot);
			Subscript *newCpivot = malloc(sizeof(Subscript));
			newCpivot->index = p_c;
			newCpivot->nonzeros = col_counts[p_c];
			newCpivot->prev = newCpivot->next = NULL;
			insertSubAtEnd(ColPivots, newCpivot);

			// 更新row_counts
			row_counts[p_r] = -1;			// use -1 to indicate the row has an elelemnt was chosen as pivot
			for (i=0; i<nrow; i++) {
				if ( (row_counts[i] != -1) && (A[i][p_c] != 0) ) {
					row_counts[i] -= 1;
				}
			}
			// 更新ColID_list
			Subscript *ss_pt = ColID_lists[col_counts[p_c]]->ssFirst;
			while (ss_pt != NULL) {
				if (ss_pt->index == p_c)
					break;
				ss_pt = ss_pt->next;
			}
			if (ss_pt->index == p_c) {
				removeSubscript(ColID_lists[col_counts[p_c]], ss_pt);
				free(ss_pt);
				ss_pt = NULL;
			}

			inactives[p_c] = 2;
			active -= 1;
			selected_pivots += 1;
			//printf("Inactivated: %d Pivots Found: %d\n", inactivated, selected_pivots);
		} else {
			// no singleton row can be found, declare one column with the most nonzeros as inactive
			// Note: other algorithms may be used to choose a column to inactivate
			for (i=max_col1s; i>=0; i--) {
				Subscript *ss_pt = ColID_lists[i]->ssFirst;
				//if ((ss_pt != NULL) && (inactives[ss_pt->index] == 0))
				if (ss_pt != NULL) {
					inactives[ss_pt->index] = 1;
					inactivated += 1;
					active -= 1;
					for (k=0; k<nrow; k++) {
						if ((A[k][ss_pt->index] != 0) && (row_counts[k] != -1))
							row_counts[k] -= 1;
					}

					removeSubscript(ColID_lists[i], ss_pt);
					free(ss_pt);
					ss_pt = NULL;
					break;
				}
			}
			selected_pivots += 1;
			//printf("Inactivated: %d Pivots Found: %d\n", inactivated, selected_pivots);
		}

	}
	// assign pivots for the dense part (any ordering is fine)
	for (i=0; i<ncolA; i++) {
		if (inactives[i] == 1) {
			// find an arbitray row that has not been selected 
			int j_candidate;
			for (j=0; j<nrow; j++) {
				if (row_counts[j] != -1) {
					j_candidate = j;
					if (A[j][i] != 0)
						break;
				}
			}
			j = j_candidate;

			// 保存该pivot的坐标
			Subscript *newRpivot = malloc(sizeof(Subscript));
			newRpivot->index = j;
			newRpivot->nonzeros = 0;
			newRpivot->prev = newRpivot->next = NULL;
			insertSubAtEnd(RowPivots, newRpivot);
			Subscript *newCpivot = malloc(sizeof(Subscript));
			newCpivot->index = i;
			newCpivot->nonzeros = 0;
			newCpivot->prev = newCpivot->next = NULL;
			insertSubAtEnd(ColPivots, newCpivot);
			row_counts[j] = -1;
			inactives[i] = 2;
			//printf("pivot at (%d, %d) is chosen for the dense part.\n", j, i);
		}
	}

	// free up memories
	free(inactives);

	for (i=0; i<max_col1s+1; i++)
		free_subscriptList(ColID_lists[i]);
	free(ColID_lists);

	free(row_counts);
	free(col_counts);
	return inactivated;
}

/*********************************************************************************************
 *      Reshape matrix
 * Given pivot sequence, re-order matrices: the i-th pivot (pr_i, pc_i) will be moved to the
 * coordiate (pr_i,pc_i) of the new matrix.
 *********************************************************************************************/
static void reshape_matrix(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[], ssList *RowPivots, ssList *ColPivots)
{
	int i, j, k;
	int row_id, col_id;
	Subscript *sub_row, *sub_col;
	Subscript *ss_pt;
	sub_row = RowPivots->ssFirst;
	sub_col = ColPivots->ssFirst;
	GF_ELEMENT temp;
	for (k=0; k<ncolA; k++) {
		row_id = sub_row->index;
		col_id = sub_col->index;

		// swap rows first
		if (row_id != k) {
			for (j=0; j<ncolA; j++) {
				temp = A[k][j];
				A[k][j] = A[row_id][j];
				A[row_id][j] = temp;
			}
			// For row operations, B needs to be swapped as well
			for (j=0; j<ncolB; j++) {
				temp = B[k][j];
				B[k][j] = B[row_id][j];
				B[row_id][j] = temp;
			}

			// due to the swap of rows, update RowPivots to track positions of the rest pivots
			ss_pt = sub_row->next;
			while (ss_pt != NULL) {
				if (ss_pt->index == k)
					ss_pt->index = row_id;

				ss_pt = ss_pt->next;
			}
		}

		// swap columns next
		if (col_id != k) {
			for (i=0; i<nrow; i++) {
				temp = A[i][k];
				A[i][k] = A[i][col_id];
				A[i][col_id] = temp;
			}
			// due to the swap of columns, update ColPivots to track positions of the rest of pivots
			ss_pt = sub_col->next;
			while (ss_pt != NULL) {
				if (ss_pt->index == k)
					ss_pt->index = col_id;

				ss_pt = ss_pt->next;
			}
		}

		sub_row = sub_row->next;
		sub_col = sub_col->next;
	}
}

// Re-order columns of the matrix
static void permute_matrix_columns(int nrow, int ncolA, GF_ELEMENT *A[], ssList *ColPivots)
{
	int i, j, k;
	int col_id;
	Subscript *sub_col;
	Subscript *ss_pt;
	sub_col = ColPivots->ssFirst;
	GF_ELEMENT temp;
	for (k=0; k<ncolA; k++) {
		col_id = sub_col->index;

		// swap columns
		if (col_id != k) {
			for (i=0; i<nrow; i++) {
				temp = A[i][k];
				A[i][k] = A[i][col_id];
				A[i][col_id] = temp;
			}
			// due to the swap of columns, update ColPivots to track positions of the rest pivots
			ss_pt = sub_col->next;
			while (ss_pt != NULL) {
				if (ss_pt->index == k)
					ss_pt->index = col_id;

				ss_pt = ss_pt->next;
			}
		}
		sub_col = sub_col->next;
	}
}

// insert a subscript object at the beginning of the list
static void insertSubAtBeginning(ssList *sub_list, Subscript *sub)
{
	if (sub_list->ssFirst == NULL) {
		sub_list->ssFirst = sub_list->ssLast = sub;
		sub->prev = NULL;
		sub->next = NULL;
	} else {
		sub_list->ssFirst->prev = sub;
		sub->next = sub_list->ssFirst;
		sub_list->ssFirst = sub;
		sub->prev = NULL;
	}
}

// insert a subscript object at the end of the list
static void insertSubAtEnd(ssList *sub_list, Subscript *sub)
{
	if (sub_list->ssFirst == NULL && sub_list->ssLast == NULL) {
		sub_list->ssFirst = sub_list->ssLast = sub;
		sub->prev = NULL;
		sub->next = NULL;
	} else {
		sub_list->ssLast->next = sub;
		sub->prev = sub_list->ssLast;
		sub->next = NULL;
		sub_list->ssLast = sub;
	}
}

// remove a subscript object from the list (no free() operation yet)
static void removeSubscript(ssList *sub_list, Subscript *sub)
{
	if (sub->prev == NULL && sub->next == NULL) {
		// only one left in the list
		sub_list->ssFirst = sub_list->ssLast = NULL;
	} else if (sub->prev == NULL && sub->next != NULL) {
		// the element is the head
		sub_list->ssFirst = sub->next;
		sub_list->ssFirst->prev = NULL;
	} else if (sub->prev != NULL && sub->next == NULL) {
		// the element is the tail
		sub_list->ssLast = sub->prev;
		sub->prev->next = NULL;
	} else {
		// the element is in the middle of the list
		sub->next->prev = sub->prev;
		sub->prev->next = sub->next;
	}
}

// free the subscript list
static void free_subscriptList(ssList *sub_list)
{
	Subscript *sub = sub_list->ssFirst;
	while (sub != NULL) {
		sub_list->ssFirst = sub->next;
		free(sub);
		sub = sub_list->ssFirst;
	}
	free(sub_list);
	sub_list = NULL;
}
