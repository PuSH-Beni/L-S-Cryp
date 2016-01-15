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


    Mat *matPlain = newMat(8, 8, 0);

    for (i = 0; i < 8; ++i) {
		int tem;
        fscanf(fin, "%x ", &tem);
		*(matPlain->vect + i) = (BYTE)tem;
    }
	
    Mat *matMul  = newMat(8, 8, 0);
    for (i = 0; i < 8; ++i) {
        fscanf(fin, "%x ", matMul->vect + i);
    }

	/* split */
	Mat **parts;
    parts = split(matPlain, 2, 1);
	for (i = 0; i < 4; ++i) {
		printf("%x ", *(parts[0]->vect + i));
	}
	printf("\n ==>Next sub-mat:\n");
	for (i = 0; i < 4; ++i) {
		printf("%x ", *(parts[1]->vect + i));
	}
    
   
	printf("\n ==>Mat:\n");	
	for (i = 0; i < 8; ++i) {
		printf("%x ", *(matPlain->vect + i));
	}

	printf("\n ==>Mat2:\n");
	for (i = 0; i < 8; ++i) {
		printf( "%x ", *(matMul->vect + i));
	}
    

	/* multiply */
	Mat *prod;
	prod = multiply(matPlain, matMul);
	printf("\n ==>MULTI_RES:\n");
	for (i = 0; i < 8; ++i) {
		printf("%02x ", *(prod->vect + i));
	}


	/* sboxes */
	/*Mat *sout;
	sout = sboxes(matPlain, parts[1]);
	printf("\n ==>Sboxes:\n");
	for (i = 0; i < 8; ++i) {
		printf("%02x ", *(sout->vect + i));
	}*/


	/* Indentity Matrix as matL */
	Mat *matI = newMat(DIM_S, DIM_L, 0);
	identMat64(matI->vect);
	printf("\n ==>Identity:\n");
	for (i = 0; i < 8; ++i) {
		printf("%02x ", *(matI->vect + i));
	}

	/* Zero Matrix */
	Mat *matZero = newMat(DIM_S, DIM_L, 0);
	Mat *matZeroHalf = newMat(DIM_S/2, DIM_L, 0);
    

    /* key */
    Mat *matK = matI;

    /* enc */

	Mat *roundOut;
	for (i = 0; i < 8; ++i){
		 roundOut = encryp(matPlain, matK, matZero, matI, matZeroHalf, ROUNDS);
		 
		 matPlain = roundOut;
	}
    

    printf("\n ==>LSout:\n");
    for (i = 0; i < 8; ++i) {
        printf("%02x ", *(matPlain->vect + i));
    }


    fclose(fin);
    fclose(fout);
    return 0;
}
