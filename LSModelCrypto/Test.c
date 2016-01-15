//
//  Test.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//	This file is a test file, to test the functions implemented in the project.


#include "Fundamentals.h"
#include "LSbox.h"
#include <stdio.h>

#define FILE_IO_TEST		1
#define CONSTRUCT_MAT_TEST  1
#define MULTIPLY_TEST		1
#define ENCRYPT_TEST		1


static FILE *fin, *fout;
static int i;

int main(){

	fin = fopen("in.txt", "r");
	fout = fopen("out.txt", "w");

	/* =============    TEST BEGINS  =================== */

#if CONSTRUCT_MAT_TEST
	Mat *matA = newMat(DIM_S, DIM_L, NULL, 0x00);
	Mat *matB = newMat(DIM_S, DIM_L, NULL, 0x00);
#endif

#if FILE_IO_TEST

	/* Dectect the input file existed or not */
	if (fin == NULL || fout == NULL) {
		printf("File Doesnt Exist\n");
		return 1;
	}
	else{
		printf("File Read Seccessfully\n");
	}

	/* Read the first matrix */
	for (i = 0; i < 8; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matA->vect + i) = (BYTE)tem;
	}
	/* Show it */
	printf("\n ==>Mat:\n");
	for (i = 0; i < 8; ++i) {
		printf("%x ", *(matA->vect + i));
	}

	/* Read the second matrix */
	for (i = 0; i < 8; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matB->vect + i) = (BYTE)tem;
	}
	/* Show it */
	printf("\n ==>Mat2:\n");
	for (i = 0; i < 8; ++i) {
		printf("%x ", *(matB->vect + i));
	}
#endif

#if MULTIPLY_TEST

	/* Matrix-Multiply */

	Mat *prod;
	prod = multiply(matA, matB);
	printf("\n ==>MULTI_RES:\n");
	for (i = 0; i < 8; ++i) {
		printf("%02x ", *(prod->vect + i));
	}
#endif

#if ENCRYPT_TEST

	/* Pre-Calculate some constant matrices */
	newPreCal();

	/* L-S-Model Eencryption */

	Mat *cipher;
	cipher = encryp(matA, matB);

	printf("\n ==>LSout:\n");
	for (i = 0; i < 8; ++i) {
		printf("%02x ", *(cipher->vect + i));
	}

#endif




	/* Free the allocated mems and Deallocate all pointers */
#if ENCRYPT_TEST
	dePostCal();
#endif

#if CONSTRUCT_MAT_TEST
	deMat(matA);
	deMat(matB);
#endif

#if FILE_IO_TEST
	fclose(fin);
	fclose(fout);
#endif


	return 0;
}
