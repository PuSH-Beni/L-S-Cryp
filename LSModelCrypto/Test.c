//
//  Test.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//	This file is a test file, to test the functions implemented in the project.



#include "LSbox.h"
#include <stdio.h>

/* <time.h> must be declared explicitly in C99 */
#include <time.h>


/*=========================================================*/
/* MARK:  Toggle Of Test Options     */
/*=========================================================*/
#define RAND_TEST						 1
#define CONSTRUCT_MAT_TEST				 1
#define FILE_IO_TEST					 1
#define TRANSPOSE_TEST					 1

#if MASK
#define RAND_ORDER_TEST				     0
#define TENSOR_PRODUCT_TEST				 0
#define GEN_RAND_MAT_TEST				 0
#define SET_UP_TEST						 0
#define ENCODE_TEST						 0
#define REFRESH_TEST					 0
#define BIT_AND_TEST					 0
#define ADD_TEST						 0
#endif

#define MULTIPLY_TEST					 1
#define ENCRYPT_TEST					 1




/*=========================================================*/
/* MARK: Main Function Begins    */
/*=========================================================*/
int main(){

	FILE *fin, *fout;
	int i, j, bts;
	Res res;
	Mat *matTem;
	fin = fopen("in.txt", "r");
	fout = fopen("out.txt", "w");

	srand((WORD)time(NULL));

	/*=====================================================*/
	/* =============    TEST BEGINS  =================== */
	/*=====================================================*/
#if RAND_TEST
      
    int origArray[] = { 0, 1, 2, 3, 4, 5, 6, 7};
    for (i = 0; i < LENGTH; ++i) {
        int j = rand() % LENGTH;
        
        int tem = origArray[j];
        origArray[j] = origArray[i];
        origArray[i] = tem;
        
    }
    for (i = 0; i < LENGTH; ++i) {
        printf("%d ",origArray[i]) ;
    }
    printf("\n");
#endif
    
    
#if CONSTRUCT_MAT_TEST
	Mat *matX = newMat(DIM_S, DIM_L, NULL, 0x00);
	Mat *matY = newMat(DIM_S, DIM_L, NULL, 0x00);
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
	for (i = 0; i < LENGTH; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matX->vect + i) = (BYTE)tem;
	}
	/* Show it */
	printf("\n ==>Mat X :\n");
	for (i = 0; i < LENGTH; ++i) {
		printf("%x ", *(matX->vect + i));
	}
	printf("\n");

	/* Read the second matrix */
	for (i = 0; i < LENGTH; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matY->vect + i) = (BYTE)tem;
	}
	/* Show it */
	printf("\n ==>Mat Y  :\n");
	for (i = 0; i < LENGTH; ++i) {
		printf("%x ", *(matY->vect + i));
	}
	printf("\n");
#endif

#if TRANSPOSE_TEST
	Mat *matTransX = transpose(matX);
	bts = bytesOfRow(matTransX->dim_col);
	BYTE *ptrOfTransX = matTransX->vect;
	printf("\n ==>Transpose_RES:\n");
	for (i = 0; i < matTransX->dim_row; ++i){
		for (j = 0; j < bts; ++j){
			printf("%02x ", *ptrOfTransX);
			++ptrOfTransX;
		}
	}
	printf("\n");

#endif

#if MULTIPLY_TEST

	/* Matrix-Multiply */

	Mat *prod;
	prod = multiply(matX, matY);
	printf("\n ==>MULTI_RES:\n");
	for (i = 0; i < LENGTH; ++i) {
		printf("%02x ", *(prod->vect + i));
	}
	printf("\n");
#endif



#if MASK

#if RAND_ORDER_TEST
	int *nums = randOrder();
	printf("\n ==>Rand_Order_RES:\n");
	for (i = 0; i < LENGTH; ++i){
		printf("%u ", nums[i]);
	}
	printf("\n");
#endif

#if GEN_RAND_MAT_TEST
	Mat **matsAs = genRandMat();
	for (i = 0; i < 3; ++i){
		printf("\n ==> GENERATE_MATS[%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matsAs[i]->vect + j));
		}
		printf("\n");
	}
	Mat *matE = multiply(matsAs[1], matsAs[2]);
	printf("\n ==> A x inv(A) :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matE->vect + j));
	}
	printf("\n");

#endif

#if TENSOR_PRODUCT_TEST
	Mat *tensorProd = tensorProduct(matX, matY);
	printf("\n ==> tenserProd  :\n", i);

	bts = bytesOfRow(tensorProd->dim_col);
	BYTE *ptrOfTenserPd = tensorProd->vect;
	for (i = 0; i < tensorProd->dim_row; ++i){
		for (j = 0; j < bts; ++j){
			printf("%02x ", *ptrOfTenserPd);
			++ptrOfTenserPd;
		}
		printf("\n");
	}
	printf("\n");
#endif

#if SET_UP_TEST

	////A = inv(A) = trans(A) = E
	//BYTE vecId1[] = INDENT_MAT;
	//BYTE vecId2[] = INDENT_MAT;
	//BYTE vecId3[] = INDENT_MAT;
	//Mat *matE1 = newMat(LENGTH, LENGTH, vecId1, 0x03);
	//Mat *matE2 = newMat(LENGTH, LENGTH, vecId2, 0x03);
	//Mat *matE3 = newMat(LENGTH, LENGTH, vecId3, 0x03);

	//matA = matE1;
	//matInvA = matE2;
	//matTransA = matE3;

	setup();

	printf("\n ==> A  :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matA->vect + j));
	}
	printf("\n");
	printf("\n ==> inv(A)  :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matInvA->vect + j));
	}
	printf("\n");
	printf("\n ==> trans(A)  :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matTransA->vect + j));
	}
	printf("\n");

	Mat *matE = multiply(matInvA, matTransA);
	printf("\n ==> A x inv(A) :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matE->vect + j));
	}
	printf("\n");

	printf("\n ==> \\hat{A}  :\n", i);
	for (j = 0; j < LENGTH ; ++j){
		for (i = 0; i < LENGTH; ++i){
			printf("%02x ", *(matHat->vect + j * LENGTH + i));
		}
		printf("\n");
	}
	printf("\n");
	printf("\n ==> \\grave{A}  :\n", i);
	for (j = 0; j < LENGTH; ++j){
		for (i = 0; i < LENGTH; ++i){
			printf("%02x ", *(matGrave->vect + j * LENGTH + i));
		}
		printf("\n");
	}
	printf("\n");
	printf("\n ==> \\acute{A}  :\n", i);
	for (j = 0; j < LENGTH; ++j){
		for (i = 0; i < LENGTH; ++i){
			printf("%02x ", *(matAcute->vect + j * LENGTH + i));
		}
		printf("\n");
	}
	printf("\n");
#endif

#if ENCODE_TEST
	/* encoded X */
	Mat **matEX =  encode(matX);
	for (i = 0; i < MASKD; ++i){
		printf("\n ==> Encoded_X [%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matEX[i]->vect + j));
		}
		printf("\n");

	}

	/* encoded Y */
	Mat **matEY = encode(matY);
	for (i = 0; i < MASKD; ++i){
		printf("\n ==> Encoded_Y [%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matEY[i]->vect + j));
		}
		printf("\n");
	}
#endif

#if REFRESH_TEST		
	/* Refreshing matEX */
	printf("\n ===> Refreshing the masks\n");
	res = refreshing(matEX);
	if(!res) printf("---> refreshed succeed\n");
	else return 1;
	for (i = 0; i < MASKD; ++i){
		printf("\n ==> Refreshed_X [%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matEX[i]->vect + j));
		}
		printf("\n");
	}

	
	Mat *matXres = decode(matEX);

	printf("\n ==> after refreshing, Represent X :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matXres->vect + j));
	}
	printf("\n");


	/* Refreshing matEY */
	printf("\n ===> Refreshing the masks\n");
	res = refreshing(matEY);
	if (!res) printf("---> refreshed succeed\n");
	else return 1;
	for (i = 0; i < MASKD; ++i){
		printf("\n ==> Refreshed_Y [%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matEY[i]->vect + j));
		}
		printf("\n");
	}
	Mat *matYres = decode(matEY);

	printf("\n ==> after refreshing, Represent Y :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matYres->vect + j));
	}
	printf("\n");
	
#endif

#if ADD_TEST
	Mat **matSum = addWithMask(matEX, matEY);
	for (i = 0; i < MASKD; ++i){
		printf("\n ==> Encoded_SUM [%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matSum[i]->vect + j));
		}
		printf("\n");

	}
	Mat *matSumRes = decode(matSum);

	printf("\n ==>  Represent Sum :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matSumRes->vect + j));
	}
	printf("\n");
#endif
#if BIT_AND_TEST
	/*
	Mat **matsX[MASKD];
	Mat **matsY[MASKD];
	for (i = 0; i < MASKD; ++i){
	matsX[i] = split(matEX[i], matX->dim_row, 1);
	matsY[i] = split(matEY[i], matY->dim_row, 1);
	}
	Mat *matRowX[] = { matsX[0][0], matsX[0][1], matsX[0][2], matsX[0][3] };
	Mat *matRowY[] = { matsY[0][0], matsY[0][1], matsY[0][2], matsY[0][3] };
	*/

	Mat **matEZ = bitAndWithMask(matEX, matEY);
	for (i = 0; i < MASKD; ++i){
		printf("\n ==> Encoded_Z [%d] :\n", i);
		for (j = 0; j < LENGTH; ++j){
			printf("%02x ", *(matEZ[i]->vect + j));
		}
		printf("\n");

	}
	Mat *matZres = decode(matEZ);

	printf("\n ==>  Represent Z :\n", i);
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(matZres->vect + j));
	}
	printf("\n");


#endif

#endif




#if ENCRYPT_TEST
#if MASK

	/* L-S-Model Eencryption */

	Mat *cipher;
	cipher = encrypto(matX, matY);

	printf("\n ==>LSout:\n");
	for (i = 0; i < LENGTH; ++i) {
		printf("%02x ", *(cipher->vect + i));
	}

#else
	/* Pre-Calculate some constant matrices */
	newPreCal();

	/* L-S-Model Eencryption */

	Mat *cipher;
	cipher = encrypto(matX, matY);

	printf("\n ==>LSout:\n");
	for (i = 0; i < LENGTH; ++i) {
		printf("%02x ", *(cipher->vect + i));
	}
	dePostCal();
#endif
#endif







	/*==============================ASAA=======================*/
	/* MARK: Free the allocated mems and Deallocate all pointers */
	/*=====================================================*/
#if CONSTRUCT_MAT_TEST
	deMat(matX);
	deMat(matY);
#endif

	return 0;
}
