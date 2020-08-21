/*
 * main.c
 *
 *  Created on: 8 במאי 2020
 *      Author: DELL
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <assert.h>
#include <time.h>
#include "spmat.h"
#include <stdbool.h>
#include <math.h>
#include <time.h>

 bool compareVector(int n, double* vector, double* nextvector){
	 int i;
	 double epsilon = 0.00001;
	 for (i=0;i<n;i++){
		 if (fabs(vector[i]-nextvector[i])>=epsilon)
			 return false;
	 }
	 return true;
 }

void createVector(int argc, double* nextVector,char *argv[],int n){
	int i,tmp;
	if (argc==4){
		srand(time(NULL));
		for (i=0;i<n;i++){
			nextVector[i] = rand();
		}
	}
	else { /*argc==5 */
		FILE *vectorFile;
		vectorFile =  fopen(argv[2], "r");
		assert(vectorFile != NULL);
		tmp = fread(nextVector, sizeof(double), n, vectorFile);
		assert(tmp==n);
		fclose(vectorFile);
	}
}

/*allocating spmat according to flag*/
spmat* allocate(char *flag, double *readrow, FILE* fileIn,int n){
	int i,tmp,j,arr[2],nnz;
	if (strcmp(flag, "-list") == 0) {
		return spmat_allocate_list(n);
	}

	/*allocating array flag==-array */
	nnz = 0;
	for (i=0;i<n;i++){ /*calculation nnz*/
		tmp = fread(readrow, sizeof(double), n, fileIn);
		assert(tmp==n);
		for(j=0;j<n;j++){
			if (readrow[j]!=0.0){
				nnz++;
			}
		}
	}
	rewind(fileIn);
	tmp = fread(arr, sizeof(int), 2, fileIn);
	assert (tmp==2);
	return spmat_allocate_array(n,nnz);
}

void writeFreeClose(FILE *fileOut,FILE *fileIn, double *vector, double *nextVector, double *readrow, spmat* A, int n){
	/*writing the matrix size, and result vector to fileOut*/
	int tmp,i=1;
	tmp = fwrite(&i,sizeof(int),1,fileOut);
	assert(tmp==1);
	tmp = fwrite(&n,sizeof(int),1,fileOut);
	assert(tmp==1);
	tmp = fwrite(nextVector,sizeof(double),n,fileOut);
	assert(tmp==n);

	/*free all allocations and closing files*/
	free(vector);
	free(nextVector);
	free(readrow);
	A->free(A);
	fclose(fileOut);
	fclose(fileIn);
}


int main(int argc, char *argv[])
{
	int arr[2], n,i,tmp, iter;
	char* flag;
	spmat *A;
	FILE *fileIn, *fileOut;
	double *vector, *nextVector, *readrow, timeSpan;
	clock_t end;
	clock_t start = clock();
	assert(argc==4 || argc==5);

	flag = argv[argc-1];
	fileIn = fopen(argv[1], "r");
	fileOut = fopen(argv[argc-2], "w");
	assert(fileIn!=NULL && fileOut!=NULL);
	tmp = fread(arr, sizeof(int), 2, fileIn);
	assert (tmp==2);
	n = arr[0];

	nextVector = (double*)calloc(n,sizeof(double));
	vector = (double*)malloc(n*sizeof(double));
	readrow = (double*)malloc(n*sizeof(double));
	assert(vector!=NULL && readrow!=NULL && nextVector!=NULL);

	createVector(argc,nextVector,argv,n);
	A = allocate(flag,readrow, fileIn, n);
	iter=0;

	/*power iteration and finding eigen vector*/
	for (i=0;i<n;i++){
		tmp = fread(readrow, sizeof(double), n, fileIn);
		assert(tmp==n);
		A->add_row(A,readrow,i);
	}

	do {
		double root=0.0, dotProduct=0.0;
		double *pointer = vector;
		vector = nextVector;
		nextVector = pointer;
		A->mult(A,vector,nextVector);
		iter++;
		for (i = 0; i<n; i++){
			dotProduct += nextVector[i]*nextVector[i];
		}
		if (dotProduct != 0){/*if dotProduct = 0 vector is 0 vector no need to normalize*/
			root = sqrt(dotProduct);
			 for (i=0;i<n;i++){
				 nextVector[i] = nextVector[i]/root;

			 }
		}
	}
	while (!compareVector(n, vector, nextVector));

	writeFreeClose(fileOut,fileIn,vector,nextVector,readrow,A,n);


	end = clock();
	timeSpan = ((double)(end-start)/CLOCKS_PER_SEC);
	printf("%s took %f\n",flag,timeSpan);
	printf("%d iterations",iter);
	return 0;
}
