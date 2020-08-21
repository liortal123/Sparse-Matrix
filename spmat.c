/*
 * spmat.c
 *
 *  Created on: 8 במאי 2020
 *      Author: DELL
 */
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "spmat.h"

typedef   double   DATA;

typedef struct linklist
{
	DATA  data;
	int col;
	struct linklist *next;
}linklist;


typedef struct _spmat_list {
	 linklist** rowlist;

} spmat_list;

typedef struct _spmat_array {
	 double* values;
	 int *colind,*rowptr, nnz;

} spmat_array;


/* Array Implementation */
/* Adds row i the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1) */
void	add_row_array(struct _spmat *A, const double *row, int row_i)
{
	int i,n,pointer;
	spmat_array *spa = (spmat_array*)A->private;
	n = A->n;
	pointer=spa->rowptr[row_i]; /*pointer symbolize where to add the next nnz to arrays colind and values */
	if (spa->nnz == pointer){
		spa->rowptr[row_i+1] = pointer;
		return;
	}
	for (i=0;i<n;i++){
		if (row[i]==0.0) {continue;}
		spa->colind[pointer] = i;
		spa->values[pointer] = row[i];
		pointer++;
	}
	spa->rowptr[row_i+1] = pointer; /*update pointer of next row */
}

/* Frees all resources used by A */
void	free_array(struct _spmat *A){
	spmat_array *spa = (spmat_array*)A->private;
	if(spa->nnz != 0){
		free(spa->colind);
		free(spa->values);
	}
	free(spa->rowptr);
	free(spa);
	free(A);
}

/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
void	mult_array(const struct _spmat *A, const double *v, double *result){
	int n,i,j,diff;
	double sum;
	spmat_array *spa = (spmat_array*)A->private;
	n = A->n;
	for (i=0;i<n;i++){
		int rowptr = spa->rowptr[i]; /*rowptr symbolize the index in the colind array */
		diff = spa->rowptr[i+1]-rowptr; /*diff symbolize how many nnz elements in row i */
		sum = 0.0;
		for (j=0;j<diff;j++){
			int col_ind = spa->colind[rowptr+j]; /*col_ind symbolize the column of the nnz in the original matrix */
			sum += v[col_ind]*spa->values[rowptr+j]; /*multipy the corresponding nnz value in row i and vector */
		}
		result[i] = sum;
	}
}


/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_array(int n,int nnz)
{
	spmat_array *spa;
	spmat *sp = (spmat*)malloc(sizeof(spmat));
	assert(sp!=NULL);
	sp->add_row = add_row_array;
	sp->free = free_array;
	sp->mult = mult_array;
	sp->n = n;

	spa = (spmat_array*)malloc(sizeof(spmat_array));
	assert(spa!=NULL);
	spa->nnz = nnz;
	spa->rowptr = (int*)calloc((n+1),sizeof(int));
	assert(spa->rowptr!=NULL);
	if (nnz != 0){
		spa->values = (double*)malloc(nnz*sizeof(double));
		spa->colind = (int*)malloc(nnz*sizeof(int));
		assert(spa->values!=NULL && spa->colind!=NULL);
	}

	sp->private = (void*)spa;
	return sp;
}



/* LIST Implementation */
/* Adds row i the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1) */
void add_row_list(struct _spmat *A, const double *row, int row_i)
{
	int j;
	linklist* tail = NULL;
	spmat_list *spl = (spmat_list*)A->private;

	for (j=0;j<A->n;j++)
	{
		if(row[j]==0.0) {
			continue;
		}
		if ((spl->rowlist[row_i])==NULL){
			spl->rowlist[row_i] =  (linklist*)malloc(sizeof(linklist));
			assert(spl->rowlist[row_i]!=NULL);
			tail = spl->rowlist[row_i];
			tail->data = row[j];
			tail->col = j;
			tail->next = NULL;
		}
		else {
			tail->next = (linklist*) malloc(sizeof(linklist));
			assert(tail->next !=NULL);
			tail->next->data = row[j];
			tail->next->col = j;
			tail =tail->next;
			tail->next = NULL;
		}


	}
}


/* Frees all resources used by A */
void free_list(struct _spmat *A)
{
	int i;
	spmat_list *spl = (spmat_list*)A->private;
	if (spl->rowlist != NULL){
		for (i=0;i<A->n;i++){
			while (spl->rowlist[i] != NULL){
				linklist *tmp = spl->rowlist[i]->next;
				free(spl->rowlist[i]);
				spl->rowlist[i] = tmp;

				}
			}
		free(spl->rowlist);
		}
	free(spl);
	free(A);

}

/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
void mult_list(const struct _spmat *A, const double *v, double *result)
{
	int i;
	double sum;
	spmat_list *spl = (spmat_list*)A->private;
	linklist* tmp;

	for (i=0;i<A->n;i++){
		sum=0.0;
		tmp = spl->rowlist[i];
		while(tmp!=NULL){
			sum += tmp->data*v[tmp->col];
			tmp = tmp->next;
		}
		result[i] = sum;
	}
}

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_list(int n)
{
	spmat_list *spl;
	spmat *sp = (spmat*)malloc(sizeof(spmat));
	assert(sp!=NULL);
	sp->add_row = add_row_list;
	sp->free = free_list;
	sp->mult = mult_list;
	sp->n = n;

	spl = (spmat_list*)malloc(sizeof(spmat_list));
	assert(spl!=NULL);
	spl->rowlist = (linklist**)calloc(n,sizeof(linklist*));
	assert(spl->rowlist!=NULL);
	sp->private = (void*)spl;
	return sp;
}

