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
#define CONSTRUCT_MAT_TEST				 1
#define FILE_IO_TEST					 1

#define TRANSPOSE_TEST					 0
#define MULTIPLY_TEST					 0
#define SPLIT_CAT_TEST					 0

#if MASK
#define TENSOR_PRODUCT_TEST				 0
#define GEN_RAND_MAT_TEST				 0
#define SET_UP_TEST						 0
#define ENCODE_TEST						 0
#define BIT_AND_TEST					 0
#define ADD_TEST						 0
#endif


#define ENCRYPT_TEST					 1

#define TIMES							 2

#if MASK
#if DIM_A
extern
Mat *matA, *matInvA, *matTransA, *matHat, *matGrave, *matAcute;
#endif
#endif
/*=========================================================*/
/* MARK: Main Function Begins    */
/*=========================================================*/
int main(){

	FILE *fin, *fout;
	int i, j, k, bts;
	Res res;
	Mat *matTem;
	fin = fopen("in.txt", "r");
	fout = fopen("out.txt", "a");

	srand((WORD)time(NULL));

	/*=====================================================*/
	/* =============    TEST BEGINS  =================== */
	/*=====================================================*/
    
#if CONSTRUCT_MAT_TEST
	Mat *matX = newMat(DIM_S, DIM_L, NULL, 0x00);
	Mat *matY = newMat(DIM_S, DIM_L, NULL, 0x00);
#endif

#if FILE_IO_TEST

	/* Dectect the input file existed or not */
	if (fin == NULL || fout == NULL) {
		printf( "File Doesnt Exist\n");
		return 1;
	}
	else{
		printf(  "File Read Seccessfully\n");
	}

	/* Read the first matrix */
	for (i = 0; i < DIM_L; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matX->vect + i) = (BASE)tem;
	}
	/* Show it */
	printf( "\n ==>Mat X :\n");
	for (i = 0; i < DIM_L; ++i) {
		printf( "%x ", *(matX->vect + i));
	}
	printf( "\n");

	/* Read the second matrix */
	for (i = 0; i < DIM_L; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matY->vect + i) = (BASE)tem;
	}
	/* Show it */
	printf( "\n ==>Mat Y  :\n");
	for (i = 0; i < DIM_L; ++i) {
		printf( "%x ", *(matY->vect + i));
	}
	printf( "\n");
#endif

#if TRANSPOSE_TEST
	Mat *matTransX = transpose(matX);
	bts = bytesOfRow(matTransX->dim_col);
	BASE *ptrOfTransX = matTransX->vect;
	printf("\n ==>Transpose_RES:\n");
	for (i = 0; i < matTransX->dim_row; ++i){
		for (j = 0; j < bts; ++j){
			printf("%02x ", *ptrOfTransX++);
		}
	}
	printf("\n");

#endif

#if MULTIPLY_TEST

	/* Matrix-Multiply */

	Mat *prod = multiply(matX, matY);
	printf("\n ==>MULTI_RES:\n");
	bts = bytesOfRow(prod->dim_col);
	BASE *ptrOfProd= prod->vect;
	for (i = 0; i < prod->dim_row; ++i){
		for (j = 0; j < bts; ++j){
			printf("%02x ", *ptrOfProd++);
		}
	}
	printf("\n");
#endif

#if SPLIT_CAT_TEST
	int slices = 2;
	Mat **splitRes = split(matX, slices, 1);
	Mat *catRes = cat(splitRes, slices, 1);
	printf("\n ==>split_RES:\n");
	BASE *ptrOfMat;
	for (k = 0; k < slices; ++k){
		ptrOfMat = splitRes[k]->vect;
		bts = bytesOfRow(splitRes[k]->dim_col);
		for (i = 0; i < splitRes[k]->dim_row; ++i){
			for (j = 0; j < bts; ++j){
				printf("%02x ", *ptrOfMat++);
			}
		}
		printf("\n");
	}
	printf("\n");
	printf("\n ==>cat_RES:\n");
	ptrOfMat = catRes->vect;
	bts = bytesOfRow(catRes->dim_col);
	for (i = 0; i < catRes->dim_row; ++i){
		for (j = 0; j < bts; ++j){
			printf("%02x ", *ptrOfMat++);
		}
	}
	printf("\n");
	
#endif


#if MASK

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
     /* MARK: TENSOR_PRODUCT */
#if DIM_A == 4
	Mat **matsX = split(matX,2,2);
	Mat **matsY = split(matY,2,2);
	Mat *catX = cat(matsX, 2, 2);
	printf("\n ==> catRESfor X :\n");
	for (j = 0; j < LENGTH; ++j){
		printf("%02x ", *(catX->vect + j));
	}
	printf("\n");
    printf("\n ==> splitRESfor X :\n");
    for (j = 0; j < LENGTH; ++j){
        printf("%02x ", *(matsX[1]->vect + j));
    }
    printf("\n");
    printf("\n ==> splitRESfor Y :\n");
    for (j = 0; j < LENGTH; ++j){
        printf("%02x ", *(matsY[0]->vect + j));
    }
    printf("\n");
	Mat *tensorProd = tensorProduct(matsX[1], matsY[1]);
#else
	Mat *tensorProd = tensorProduct(matX, matY);
#endif
	printf("\n ==> tenserProd  :\n");

	bts = bytesOfRow(tensorProd->dim_col);
	BASE *ptrOfTenserPd = tensorProd->vect;
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
    /* MARK: SET_UP */
	setup();

	printf("\n ==> A  :\n");
	for (j = 0; j < DIM_A; ++j){
		printf("%02x ", *(matA->vect + j));
	}
	printf("\n");
	printf("\n ==> inv(A)  :\n");
	for (j = 0; j < DIM_A; ++j){
		printf("%02x ", *(matInvA->vect + j));
	}
	printf("\n");
	printf("\n ==> trans(A)  :\n");
	for (j = 0; j < DIM_A; ++j){
		printf("%02x ", *(matTransA->vect + j));
	}
	printf("\n");

	matE = multiply(matInvA, matTransA);
	printf("\n ==> A x inv(A) :\n");
	for (j = 0; j < DIM_A; ++j){
		printf("%02x ", *(matE->vect + j));
	}
	printf("\n");

	printf("\n ==> \\hat{A}  :\n");
    bts = bytesOfRow(matHat->dim_col);
    BASE *ptrOfHat = matHat->vect;
    for (i = 0; i < matHat->dim_row; ++i){
        for (j = 0; j < bts; ++j){
            printf("%02x ", *ptrOfHat);
            ++ptrOfHat;
        }
        printf("\n");
    }
    printf("\n");
    
	printf("\n");
	printf("\n ==> \\grave{A}  :\n");
    bts = bytesOfRow(matGrave->dim_col);
    BASE *ptrOfGrave = matGrave->vect;
    for (i = 0; i < matGrave->dim_row; ++i){
        for (j = 0; j < bts; ++j){
            printf("%02x ", *ptrOfGrave);
            ++ptrOfGrave;
        }
        printf("\n");
    }
    printf("\n");
	printf("\n ==> \\acute{A}  :\n");
    bts = bytesOfRow(matAcute->dim_col);
    BASE *ptrOfAcute = matAcute->vect;
    for (i = 0; i < matAcute->dim_row; ++i){
        for (j = 0; j < bts; ++j){
            printf("%02x ", *ptrOfAcute);
            ++ptrOfAcute;
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
     /* MARK: BITAND */
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

#endif /* MASK */




#if ENCRYPT_TEST
     /* MARK: ENCRYPTO */

	/* L-S-Model Eencryption */

	Mat *cipher;
	printf( "\n ==> begins\n");

#if DIM_A
	setup();
#endif
	newPreCal();
	double time_Start = (double)clock();

	for (j = 0; j != TIMES; ++j){
		cipher = encrypto(matX, matY);

		if (j < 2){
			printf("\n==>LSout:\n");
			for (i = 0; i != DIM_L; ++i) {
				printf("%02x ", *(cipher->vect + i));
			}
		}
		deMat(cipher);
		
	}

	double time_End = (double)clock();
	dePostCal();
	printf("\n ===> done\n");
	printf("%.fms\n", (time_End - time_Start));
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
