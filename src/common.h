/*--------------------- common.h ------------------------
 *  Internal header file 
 *------------------------------------------------------*/
#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#ifndef TRACE_LEVEL 
#define TRACE_LEVEL	4
#endif
/* log levels */
#define LOG_ERROR		1
#define LOG_WARNING		2
#define LOG_INFO		3
#define LOG_DEBUG		4
#define LOG_TRACE		5

#define ALIGN(a, b) ((a) % (b) == 0 ? (a)/(b) : (a)/(b) + 1)
#define RESIDUAL(a, b) ((b) * ALIGN((a), (b)) - (a))

// node of singly linked list
struct node {
    int data;								
    struct node *next;					
};

struct node_list {
    struct node *first;
    struct node *last;
};

/* common.c */
int has_item(int array[], int item, int length);
void append_to_list(struct node_list *list, struct node *nd);
int remove_from_list(struct node_list *list, int data);
int exist_in_list(struct node_list *list, int data);
void clear_list(struct node_list *list);
void free_list(struct node_list *list);
#endif /* COMMON_H */
