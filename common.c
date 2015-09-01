#include "common.h"



// check if an item is existed in an int array
int has_item(int array[], int item, int length)
{
	int index = -1;
	for (int i=0; i<length; i++) {
		if (item == array[i]) {
			index = i;
			break;
		}
	}
	return index;
}

void append_to_list(struct node_list *list, struct node *nd)
{
	if (list->first == NULL)
		list->first = list->last = nd;
	else {
		list->last->next = nd;
		list->last = nd;
	}
}

// Remove the first node whose data is equal to "data"
// Note: this function should only be used in applications
//       where nodes in the list have unique data
int remove_from_list(struct node_list *list, int data)
{
	struct node *prev = NULL;
	struct node *curr = list->first;
	while (curr != NULL) {
		if (curr->data == data) {
			// shorten list
			if (curr == list->first && curr == list->last) { 					// list contains only one node
				list->first = list->last = NULL; 					
			} else if (curr == list->first && curr != list->last) {				// head node is to be removed
				list->first = curr->next;
			} else if (curr != list->first && curr == list->last) {				// tail node is to be removed
				list->last = prev;
				list->last->next = NULL;
			} else {
				prev->next = curr->next;
			}

			free(curr);
			curr = NULL;
			return 0;
		}
		prev = curr;
		curr = curr->next;
	}
	return -1;
}

int exist_in_list(struct node_list *list, int data)
{
	struct node *p = list->first;
	while (p != NULL) {
		if (p->data == data)
			return 1;
		p = p->next;
	}
	return 0;
}

// clear nodes in a list, but keep the list structure alive
void clear_list(struct node_list *list)
{
	struct node *nd = list->first;
	struct node *ndNext;
	while (nd != NULL) {
		ndNext = nd->next;
		free(nd);
		nd = ndNext;
	}
	list->first = list->last = NULL;
}

// Free a list, which include clear nodes in a list and free
// the list structure in the end.
void free_list(struct node_list *list)
{
	clear_list(list);
	free(list);
}
