//
//  Test.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//
#include "Fundamentals.h"
#include "LSbox.h"
#include <stdio.h>

static FILE *fin, *fout;

int main(){


    fin = fopen("in.txt","r");
    fout = fopen("out.txt","w");

    if (fin == NULL || fout == NULL) {
        printf("File Doesnt Exist\n");
        return 1;
    }

    printf("File Read Seccessfully\n");

    int i;
	newPreCal();

    Mat *matA = newMat(DIM_S, DIM_L, NULL, 0x00);
	Mat *matB = newMat(DIM_S, DIM_L, NULL, 0x00);

    for (i = 0; i < 8; ++i) {
		int tem;
        fscanf(fin, "%x ", &tem);
		*(matA->vect + i) = (BYTE)tem;
    }

	for (i = 0; i < 8; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matB->vect + i) = (BYTE)tem;
	}
	   
	printf("\n ==>Mat:\n");	
	for (i = 0; i < 8; ++i) {
		printf("%x ", *(matA->vect + i));
	}

	printf("\n ==>Mat2:\n");
	for (i = 0; i < 8; ++i) {
		printf( "%x ", *(matB->vect + i));
	}
    

	/* multiply */
	Mat *prod;
	prod = multiply(matA, matB);
	printf("\n ==>MULTI_RES:\n");
	for (i = 0; i < 8; ++i) {
		printf("%02x ", *(prod->vect + i));
	}

	


	/* LS ENC */

	Mat *cipher;
	cipher = encryp(matA, matB);
    printf("\n ==>LSout:\n");
    for (i = 0; i < 8; ++i) {
        printf("%02x ", *(cipher->vect + i));
    }

	dePostCal();
    fclose(fin);
    fclose(fout);
    return 0;
}
