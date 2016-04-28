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
extern
BYTE matA[DIM_A], matInvA[DIM_A], matTransA[DIM_A];
#endif
/*=========================================================*/
/* MARK: Main Function Begins    */
/*=========================================================*/
int main(){
	Res res;
	const BYTE plainT[DIM_L]  = { 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf };
	const BYTE cipherK[DIM_L] = { 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf };
	const int dims[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	BYTE cipher[DIM_L] = { 0 };
	srand((unsigned)time(NULL));
#if DIM_A
	setupEnc();
	getMatT();
#endif

	BYTE x[8] = {0xc1, 0x4a, 0xe0, 0xd8, 0xc9, 0x8e, 0x69, 0xeb};
	BYTE y[8] = {0xe1, 0x2a, 0xeb, 0x53, 0x88, 0x8f, 0xe2, 0x4b};
	int dimsA[4] = { 8, 8, 8, 8 };
	BYTE pd[8] = { 0 };
	BYTE trans[8] = { 0 };
	res = transpose(trans, y, dimsA);
	res = multiply(pd, matInvA, matTransA, dimsA);
	

	double timeStart = (double)clock();
	res = encrypto(cipher, plainT, cipherK);
	double timeEnd = (double)clock();
	printf("%.fms",timeEnd - timeStart);
	return 0;
}
	
