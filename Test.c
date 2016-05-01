//
//  Test.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//	This file is a test file, to test the functions implemented in the project.



#include "LSbox.h"
#include <stdio.h>
#include <time.h>


/*=========================================================*/
/* MARK:  Toggle Of Test Options     */
/*=========================================================*/

#define ENCRYPT_TEST					 1
#define TIMES							 2
#if DIM_A
#if DIM_A
extern
BYTE matA[DIM_A], matInvA[DIM_A], matTransA[DIM_A];
#endif

#if !DIVIDE && DIM_L == 16
extern
BYTE matAs[L_SIZE], matInvAs[L_SIZE], matTransAs[L_SIZE];
#endif
#endif
/*=========================================================*/
/* MARK: Main Function Begins    */
/*=========================================================*/
void outputMat(const BYTE *mat, int byts){
	int i;
	for (i = 0; i < byts; ++i){
		printf("%02x ", mat[i]);
	}
}
int main(){
	Res res;
	const BYTE plainT[DIM_L]  = { 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf };
	const BYTE cipherK[DIM_L] = { 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf };
	const int dims[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	BYTE cipher[DIM_L] = { 0 };
	srand((unsigned)time(NULL));

#if DIM_A
	res = setupEnc();
	CHECK(res);
	res = getMatT();
	CHECK(res);
#endif

	BYTE x[8] = {0xc1, 0x4a, 0xe0, 0xd8, 0xc9, 0x8e, 0x69, 0xeb};
	BYTE y[8] = {0xe1, 0x2a, 0xeb, 0x53, 0x88, 0x8f, 0xe2, 0x4b};
	BYTE E[8] = UNIT_MAT;
	int dimsA[4] = { 8, 8, 8, 8 };
	BYTE pd[16*5] = { 0 };
	BYTE encoded[16*5] = { 0 };
	BYTE trans[8 * 2] = { 0 };
	
	printf("  --> PlainText:\n");
	outputMat(plainT, DIM_L);
	printf("\n\n  --> CipherKey:\n");
	outputMat(cipherK, DIM_L);

#if DIM_A
	res = transpose(trans, matA, dimsA);
	//test the matrix A
	printf("\n\n  --> matA:\n");
	outputMat(matA, 8);
	printf("\n\n  --> matInvA:\n");
	outputMat(matInvA, 8);
	printf("\n\n  --> matA x matInvA:\n");
	res = multiply(pd, matInvA, matTransA, dimsA);
	outputMat(pd, 8);

	// test the matrix Big_A
	printf("\n\n  --> matAs:\n");
	outputMat(matAs, 32);
	printf("\n\n  --> matInvAs:\n");
	outputMat(matInvAs, 32);
	int dimsL[4] = {16,16,16,16};
	printf("\n\n  --> matAs x matInvAs:\n");
	res = multiply(pd, matInvAs, matTransAs, dimsL);
	outputMat(pd, 32);
#endif
#if MASK
	// test the bitAndWithMask
	res = encode(encoded, plainT);
	printf("\n\n  --> Encoded PlainText:\n");
	outputMat(encoded, 8 * 4);

	res = decode(trans,encoded);
	printf("\n\n  --> Decoded PlainText:\n");
	outputMat(trans, 16);

	int dimsBA[4] = { 1, DIM_L, 1, DIM_L };
	res = bitAnd(trans, plainT, plainT+2, dimsBA);
	printf("\n\n  --> Original BitAnd:\n");
	outputMat(trans, 2);

	res = bitAndWithMask(pd, encoded, encoded + 2, dimsBA);
	printf("\n\n  --> BitAnd With Mask:\n");
	outputMat(trans, 2);
#endif
	// test the Encrypto
	double timeStart = (double)clock();
	res = encrypto(cipher, plainT, cipherK);
	res = encrypto_fixed();
	double timeEnd = (double)clock();
	printf("\n\n  --> result:\n");
	outputMat(cipher, 16);

	// output the time 
	printf("\n\n  --> timeCost:\n");
	printf("  %.fms",timeEnd - timeStart);

	return 0;
}
	
