//
//  Fundamentals.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//


#include "Fundamentals.h"
 

/*=========================================================*/
/*   MARK: Declarition     */
/*=========================================================*/
#if DIM_A
/* DIM_A == 4 or DIM_A == 8 */
Mat *matA, *matInvA, *matTransA, *matHat, *matGrave, *matAcute;
#if !DIVIDE
Mat *matAs, *matInvAs, *matTransAs;
#endif

#endif

/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/


void X_aligned_memcpy_sse2(void* dest, const void* src, const unsigned long lengths)
{

	__asm
	{
		mov esi, src;    //src pointer
		mov edi, dest;   //dest pointer
		mov ebx, lengths;   //ebx is our counter 
		shr ebx, 7;      //divide by 128 (8 * 128bit registers)


	loop_copy:
		prefetchnta 128[ESI]; //SSE2 prefetch
		prefetchnta 160[ESI];
		prefetchnta 192[ESI];
		prefetchnta 224[ESI];

		movdqa xmm0, 0[ESI]; //move data from src to registers
		movdqa xmm1, 16[ESI];
		movdqa xmm2, 32[ESI];
		movdqa xmm3, 48[ESI];
		movdqa xmm4, 64[ESI];
		movdqa xmm5, 80[ESI];
		movdqa xmm6, 96[ESI];
		movdqa xmm7, 112[ESI];

		movntdq 0[EDI], xmm0; //move data from registers to dest
		movntdq 16[EDI], xmm1;
		movntdq 32[EDI], xmm2;
		movntdq 48[EDI], xmm3;
		movntdq 64[EDI], xmm4;
		movntdq 80[EDI], xmm5;
		movntdq 96[EDI], xmm6;
		movntdq 112[EDI], xmm7;

		add esi, 128;
		add edi, 128;
		dec ebx;

		jnz loop_copy; //loop please
	loop_copy_end:
	}
}

/* Caculate the num of bytes in one row */
#if !TEST
static
#endif
int bytesOfRow(
        int col
        )
{
    // bytes of each row :: if dim_col < LENGTH, allocate a byte as well
    if (col < LENGTH) return 1;
    int bytes;
#if LENG16
    bytes = col >> 4;
#elif LENG8
    bytes = col >> 3; // col / LENGTH
#endif
    bytes = (col % LENGTH ) ? bytes + 1 : bytes;
    return bytes;
}




/* Shift bit form j -> i */
#if !TEST
static
#endif
BASE shiftBit(
        BASE orig,
        int i,
        int j
        )
{
    BASE tem = (IDENT >> j) & orig;

    if (j > i) tem <<= (j - i);
    else if (i > j) tem >>= (i - j);
    else return tem;

    return tem;
}



/* Get sum of a byte */

#if !TEST
static
#endif
BASE sumOfByte(
        BASE byte
        )
{
    /* Using the fast method to get the numbers of '1' in a byte */
    BASE ret = 0;
    while (byte > 0){
        ret ^= IDENT;
        byte &= (byte - 1);
    }

    return ret;
}




/* Generate random order */
#if !TEST
static
#endif
int *randOrder(
        )
{
    /* Initialize random seed */
    /* Just do this in the main() func */
    //srand((WORD)time(NULL));

    int *nums = (int *)malloc(LENGTH * sizeof(int));
    if (nums == NULL)return NULL;
    int i;

    for (i = 0; i < LENGTH; ++i){
        nums[i] = i;
    }

    for (i = 0; i < LENGTH; ++i){
        int j = rand() % LENGTH;
        int tem = nums[j];
        nums[j] = nums[i];
        nums[i] = tem;
    }
    return nums;
}



/* Transpositon */
#if !TEST
static
#endif
Mat *transpose(
        const Mat *matOrig
        )
{
    if (matOrig == NULL) return NULL;

    /* Memory allocated */
    Mat *matRet = newMat(matOrig->dim_col, matOrig->dim_row, NULL, 0x00);
    if (matRet == NULL) return NULL;

    int bytesOfRowOrig = bytesOfRow(matOrig->dim_col);
    int bytesOfRowRet = bytesOfRow(matOrig->dim_row);

    /* Transposing */
    int colOrig, rowOrig;
    int cntBytesOrig, cntBytesRet;
    int offsetOrig, offsetRet;
    BASE vectTransed;
    BASE *ptrOfVectRet, *ptrOfVectOrig;
    for (colOrig = 0; colOrig < matOrig->dim_col; ++colOrig){
        cntBytesOrig = colOrig / LENGTH;
        offsetOrig = colOrig % LENGTH;
        for (rowOrig = 0; rowOrig < matOrig->dim_row; ++rowOrig){
            /* the bit (rowOrig, colOrig) */
            cntBytesRet = rowOrig / LENGTH;
            offsetRet = rowOrig % LENGTH;

            ptrOfVectOrig = matOrig->vect + bytesOfRowOrig * rowOrig + cntBytesOrig;
            ptrOfVectRet = matRet->vect + bytesOfRowRet * colOrig + cntBytesRet;

            vectTransed = shiftBit(*ptrOfVectOrig, offsetRet, offsetOrig);
            (*ptrOfVectRet) ^= vectTransed;
        }
    }
    return matRet;

}





/* ============================================================== */
/*      MARK: Functions For Masked Model                  */
/* ============================================================== */
#if MASK


/* Generate DIM_A x DIM_A random matrices A, inv(A), trans(A)  */
#if DIM_A
#if !TEST
static
#endif
Mat **genRandMat(
        )
{

    /* Memory allocated */
    Mat **retMats = (Mat **)malloc(3 * sizeof(Mat *));
    if (retMats == NULL) return NULL;
    BASE vecId[] = IDENT_MAT;
    Mat *matE = newMat(DIM_A, DIM_A, vecId, 0x03);


    Mat *matATem = newMat(DIM_A, DIM_A, NULL, 0x00);
    Mat  *matTransP = NULL, *matTem;


    /*  get matrix A and inv(A) temporarily */
    // int *randRow = randOrder();
    int i;
    for (i = 0; i < DIM_A; ++i) matATem->vect[i] = matE->vect[i];
    Mat *matInvATem = transpose(matATem);
    Mat *matTransATem = transpose(matATem);
    deMat(matATem);

    /* get matrix P */
    //srand((WORD)time(NULL));
    for (i = 0; i < DIM_A; ++i){

        BASE rowToAdd = (BASE)rand();

        BASE zeroToSet = IDENT >> i;
        rowToAdd &= (~zeroToSet);
#if DIM_A == 4
        rowToAdd &= 0xf0;
#endif
        matE->vect[i] ^= rowToAdd;

        /* Tanspose(P) */
        matTem = matTransP;
        matTransP = transpose(matE);
        deMat(matTem);

        /*  A = P x A   ==>  A^T = A^T x P^T  */
        matTem = matTransATem;
        matTransATem = multiply(matTransATem, matE);
        deMat(matTem);

        /*  A^{-1} = A^{-1} x P     */
        matTem = matInvATem;
        matInvATem = multiply(matInvATem, matTransP);
        deMat(matTem);

        /* Refresh matE to Identity Matrix */
        matE->vect[i] ^= rowToAdd;
    }

    retMats[2] = matTransATem;
    retMats[1] = matInvATem;
    retMats[0] = transpose(matTransATem);


    /* deallocate */
    //free(randRow);
    deMat(matTransP);
    deMat(matE);

    return retMats;
}




/* Generate hatA, graveA and acuteA */
/* All of them are DIM_A x DIM_A^2 */
#if !TEST
static
#endif
Mat *hatA(
        )
{
    if (matTransA == NULL || matInvA == NULL) return NULL;

    Mat *matRight = newMat(DIM_A * DIM_A, DIM_A, NULL, 0x00);
    if (matRight == NULL) return NULL;

    /* get transposition of the right matrix: (E' x (A tp B))^T */
    int i, j;
    BASE *ptrOfMatATi;
    BASE *ptrOfMatATj, *ptrOfMatR;
    for (i = 0; i < DIM_A; ++i){
        ptrOfMatATi = matTransA->vect + i;
        for (j = 0; j < DIM_A; ++j){
            ptrOfMatATj = matTransA->vect + j;
            ptrOfMatR = matRight->vect + (i * DIM_A + j);
            (*ptrOfMatR) = (*ptrOfMatATi) & (*ptrOfMatATj);
        }
    }

    /* get \hat{A} now */
    Mat *matRet = multiply(matInvA, matRight);
    deMat(matRight);

    return matRet;
}



#if !TEST
static
#endif
Mat *acuteA(
        )
{
    if (matTransA == NULL || matInvA == NULL) return NULL;

    Mat *matRight = newMat(DIM_A * DIM_A, DIM_A, NULL, 0x00);
    if (matRight == NULL) return NULL;



    int i, j;
    BASE *ptrOfMatATj, *ptrOfMatR, *ptrOfMatI;
    BASE vecId[] = IDENT_MAT;
    Mat *matIdent = newMat(DIM_A, DIM_A, vecId, 0x03);
    if (matIdent == NULL) return NULL;


    /* get transposition of the right matrix: (E' x (E tp A))^T */
    for (i = 0; i < DIM_A; ++i){
        ptrOfMatI = matIdent->vect + i;
        for (j = 0; j < DIM_A; ++j){
            ptrOfMatATj = matTransA->vect + j;
            ptrOfMatR = matRight->vect + (i * DIM_A + j);
            (*ptrOfMatR) = (*ptrOfMatI) & (*ptrOfMatATj);
        }
    }

    /* get \grave{A} now */
    Mat *matRet = multiply(matInvA, matRight);

    deMat(matRight);
    deMat(matIdent);
    return matRet;
}



#if !TEST
static
#endif
Mat *graveA(
        )
{

    if (matTransA == NULL || matInvA == NULL) return NULL;

    Mat *matRight = newMat(DIM_A * DIM_A, DIM_A, NULL, 0x00);
    if (matRight == NULL) return NULL;


    int i, j;
    BASE *ptrOfMatATi, *ptrOfMatR, *ptrOfMatI;
    /* Identity matrix */
    BASE vecId[] = IDENT_MAT;
    Mat *matIdent = newMat(DIM_A, DIM_A, vecId, 0x03);
    if (matIdent == NULL) return NULL;


    /* get transposition of the right matrix: (E' x (A tp E))^T */
    for (i = 0; i < DIM_A; ++i){
        ptrOfMatATi = matTransA->vect + i;
        for (j = 0; j < DIM_A; ++j){
            ptrOfMatI = matIdent->vect + j;
            ptrOfMatR = matRight->vect + i * DIM_A + j;
            (*ptrOfMatR) = (*ptrOfMatI) & (*ptrOfMatATi);
        }
    }

    /* get \Acute{A} now */
    Mat *matRet = multiply(matInvA, matRight);

    deMat(matRight);
    deMat(matIdent);
    return matRet;
}

#endif /* DIM_A != 0 */


/* Tensor Product(kron) for two vectors */
/* Return a n x n^2 matrix */
#if !TEST
static
#endif
Mat *tensorProduct(
        const Mat *matX,
        const Mat *matY
        )
{
    if (matX == NULL || matY == NULL) return NULL;
    if ((matX->dim_row != matY->dim_row) ||
            (matX->dim_col != matY->dim_col) ||
            (matX->dim_col > 8) )return NULL;

    int i, j;
    Mat *matR = newMat(matX->dim_row, matX->dim_col * matX->dim_col, NULL, 0x00);
    if (matR == NULL) return NULL;

#if DIM_A == 4
    int bts = bytesOfRow(matR->dim_col);
    for (i = 0; i < matX->dim_row; ++i) {
        BASE ident = IDENT;

        for (j = 0; j < matX->dim_col ; j += 2){
            BASE tem = 0x00;
            if ((BYTE)matX->vect[i] & ident){
                tem ^= matY->vect[i];
            }
            ident = ident >> 1;
            if ((BYTE)matX->vect[i] & ident){
                tem ^= (BYTE)matY->vect[i] >> 4;
            }
            matR->vect[i * bts + j / 2] = tem;
            ident >>= 1;

        }
    }
#else /* DIM_A == 8 or 0 */

    for (i = 0; i != matX->dim_row; ++i){
        BASE ident = IDENT;
        for (j = 0; j != matY->dim_col; ++j){
            matR->vect[i * matX->dim_col + j] = (matX->vect[i] & ident) ? matY->vect[i] : 0x00;
            ident >>= 1;
        }
    }


    //int btsOfMat = matX->dim_row * bytesOfRow(matX->dim_col);
    //int btsOfRowRet = bytesOfRow(matR->dim_col);
    //int btsOfRowY = bytesOfRow(matY->dim_col);
    //   for (i = 0; i != matX->dim_row; ++i){
    //       BASE ident = IDENT;
    //	int offsetOfRet = i * btsOfRowRet;
    //       for (j = 0; j != matY->dim_col; ++j){
    //		int k;
    //		for (k = 0; k != btsOfRowY; ++k){
    //			matR->vect[offsetOfRet + j * btsOfRowY + k] = (matX->vect[i] & ident) ? matY->vect[k] : 0x00;
    //		}
    //
    //           ident >>= 1;
    //       }
    //   }

    /*for (i = 0; i != matX->dim_row; ++i){
      BYTE pos = IDENT;
      for (j = 0; j != matX->dim_col; ++j){
      int btsOfRow = bytesOfRow(matX->dim_col);
      int indexOfBt = j >> LENGTH;
      pos >>= (j % LENGTH);
      if (matX->vect[indexOfBt] & pos) memcpy(matR->vect + (j * btsOfRow), matY->vect, btsOfRow);
      else memset(matR->vect + (j * btsOfRow), 0, btsOfRow);
      }
      }
      */
#endif /* DIM_A */
    return matR;
}

#endif /* MASK */


/*=========================================================*/
/*   MARK: Public   Functions      */
/*=========================================================*/

/* Construct a new Mat instance */
Mat *newMat(
        int dim_row,
        int dim_col,
        BASE *addr,
        BASE flags /*   flags :
                    *   '0x00': norm
                    *   '0x01': unit matrix or indentit matrix
                    *   '0x02': error
                    *   '0x03': 'vect' points to an existing array
                    */
        )
{
    if (dim_row <= 0 || dim_col <= 0)return NULL;

    Mat *retMat = (Mat *)malloc(sizeof(Mat));
    if (retMat == NULL) return NULL;

    /* Memory allocation for vector-array */
    int bytesOfMat = dim_row * bytesOfRow(dim_col);
    if (addr == NULL){
        retMat->vect = (BASE *)malloc(bytesOfMat * sizeof(BASE));
        if (retMat->vect == NULL)
        {
            free(retMat);
            return NULL;
        }
        memset(retMat->vect, 0, bytesOfMat * sizeof(BASE));
        retMat->flags = flags;
    }
    else{
        retMat->vect = addr;
        retMat->flags = 0x03;
    }



    /* Initialize */
    retMat->dim_row = dim_row;
    retMat->dim_col = dim_col;


    return retMat;
}


/* Deconstruct a Mat instance */
void deMat(
        Mat *matrix
        )
{
    if (matrix != NULL)
    {
        if (matrix->flags != 0x03){
            free(matrix->vect);
            matrix->vect = NULL;
        }
        free(matrix);
        matrix = NULL;
    }
}


/* Deconstruct n Mat instances */
void deMats(
        Mat **matrices,
        int cnt
        )
{
    if (matrices != NULL)
    {
        int i;
        for (i = 0; i < cnt; ++i){
            if (matrices[i] != NULL){
                if (matrices[i]->flags != 0x03){
                    free(matrices[i]->vect);
                    matrices[i]->vect = NULL;
                }
                free(matrices[i]);
                matrices[i] = NULL;
            }
        }
        matrices = NULL;
    }
}




/* ============================================================== */
/*      Functions For Masked Model                  */
/* ============================================================== */

#if MASK

#if DIM_A
/* Set-up */
void setup(
        )
{
    /* Get Matrix A, inv(A), trans(A) */
    Mat **mats = genRandMat();

    matA = mats[0];
    matInvA = mats[1];
    matTransA = mats[2];

    /* Get Matrix \hat{A}, \grave{A}, \acute{A} */
    matHat = hatA();
    matGrave = graveA();
    matAcute = acuteA();
#if !DIVIDE
    /*| inv(A), 0 |
      | 0, inv(A) | */

    matAs = newMat(DIM_L, DIM_L, NULL, 0x00);
    matInvAs = newMat(DIM_L, DIM_L, NULL, 0x00);
    matTransAs = newMat(DIM_L, DIM_L, NULL, 0x00);

    int btsOfRow = bytesOfRow(matAs->dim_col);
    int btsOfMat = matA->dim_row * bytesOfRow(matA->dim_col);
    int indexOfSlice, indexOfByte, indexOfDest = 0;
    for (indexOfSlice = 0; indexOfSlice != DIVIDE_PARTS; ++indexOfSlice){
#if DIM_A == 4
        BYTE oddFlag = (BYTE)indexOfSlice & 0x01;
#endif
        for (indexOfByte = 0; indexOfByte != btsOfMat; ++indexOfByte){
#if DIM_A == 4
            matAs->vect[indexOfDest]      = oddFlag ? matA->vect[indexOfByte] >> 4
                : matA->vect[indexOfByte];
            matInvAs->vect[indexOfDest]   = oddFlag ? matInvA->vect[indexOfByte] >> 4
                : matInvA->vect[indexOfByte];
            matTransAs->vect[indexOfDest] = oddFlag ? matTransA->vect[indexOfByte] >> 4
                : matTransA->vect[indexOfByte];
#else
            matAs->vect[indexOfDest]	  = matA->vect[indexOfByte];
            matInvAs->vect[indexOfDest]	  = matInvA->vect[indexOfByte];
            matTransAs->vect[indexOfDest] = matTransA->vect[indexOfByte];
#endif /* DIM_A == 4 */
            indexOfDest += btsOfRow;
        }
#if DIM_A == 4
        indexOfDest += oddFlag ? 1 : 0;
#else
        indexOfDest += 1;
#endif /* DIM_A == 4*/
    }


#endif /* !DIVIDE */
}

#endif /* DIM_A != 0 */

/* Encode(Mask) the plain text */
Mat **encode(
        const Mat *matPlain
        )
{
    if ( matPlain == NULL) return NULL;

    Mat **matsRet = (Mat **)malloc(MASKD * sizeof(Mat *));
    if (matsRet == NULL) return NULL;

    int i, j;
    int btsOfRow = bytesOfRow(matPlain->dim_col);
    int btsOfMat = matPlain->dim_row * btsOfRow;
    BASE *randMat = (BASE *)malloc(btsOfMat  * sizeof(BASE));
    if (randMat == NULL) return NULL;

    Mat *matSumOfRest = NULL;
    Mat *matTem = NULL;

    for (i = 1; i < MASKD; ++i){

        matsRet[i] = newMat(matPlain->dim_row, matPlain->dim_col, NULL, 0x00);
        if (matsRet[i] == NULL) return NULL;

        for (j = 0; j < btsOfMat; ++j){
            randMat[j] = (BASE)rand();
            if(matPlain->dim_col == 4) randMat[j] &= 0xf0;
        }
        memcpy(matsRet[i]->vect, randMat, btsOfMat * sizeof(BASE));

        if (matsRet[i] == NULL) return NULL;

        if (i == 1) matSumOfRest = matsRet[i];
        else {
            matTem = matSumOfRest;
            matSumOfRest = add(matSumOfRest, matsRet[i]);
            if(i > 2) deMat(matTem);
        }

    }

#if DIM_A

#if DIM_L == 16 && !DIVIDE
    matTem = add(matSumOfRest, matPlain);
    /*  = (x^T) x inv(A) ^ T */
    matsRet[0] = multiply(matTem, matInvAs);
    deMat(matTem);
#else
    matTem = add(matSumOfRest, matPlain);
    /*  = (x^T) x inv(A) ^ T */
    matsRet[0] = multiply(matTem, matInvA);
    deMat(matTem);
#endif

#else /* DIM_A == 0 */
    matsRet[0] = add(matSumOfRest, matPlain);
#endif
    if(MASKD > 2) deMat(matSumOfRest);
    free(randMat);
    return matsRet;
}


/* Decode(Unmask) the Secret text */
Mat *decode(
        Mat **matsSecret
        )
{
    if (*matsSecret == NULL) return NULL;

    int i;
    Mat *matOrig;
#if DIM_A

#if !DIVIDE && DIM_L == 16
    Mat *matSum = multiply(matsSecret[0], matAs);
    deMat(matsSecret[0]);
#else
    Mat *matSum = multiply(matsSecret[0], matA);
    deMat(matsSecret[0]);
#endif

#else /* DIM_A == 0 */
    Mat *matSum = matsSecret[0];
#endif
    for (i = 1; i < MASKD; ++i){
        matOrig = matSum;
        matSum = add(matSum, matsSecret[i]);
        deMat(matOrig);
        deMat(matsSecret[i]);
    }

    return matSum;
}


///* Refresh the old mask with a new mask */
//Res refreshing(
//               Mat **matsSecret
//               )
//{
//    if (*matsSecret == NULL) return RES_INVALID_POINTER;
//    Mat **matsRef = (Mat **)malloc(MASKD * sizeof(Mat *));
//    if (matsRef == NULL) return RES_INVALID_POINTER;
//
//    int i, j;
//    BASE *randMat = (BASE *)malloc((*matsSecret)->dim_row * sizeof(BASE));
//    if (randMat == NULL) return RES_INVALID_POINTER;
//
//    /* a ZERO matrix */
//    BASE zeroV[] = ZERO_MAT;
//
//    Mat *matSumOfRest = newMat((*matsSecret)->dim_row, (*matsSecret)->dim_col, zeroV, 0x03);
//    if (matSumOfRest == NULL) return RES_INVALID_POINTER;
//
//    Mat *matTem = NULL;
//
//    for (i = 1; i < MASKD; ++i){
//        matsRef[i] = newMat((*matsSecret)->dim_row, (*matsSecret)->dim_col, NULL, 0x00);
//        if (matsRef[i] == NULL) return RES_INVALID_POINTER;
//
//        for (j = 0; j < (*matsSecret)->dim_row; ++j){
//            randMat[j] = (BASE)rand();
//#if DIM_A == 4
//        randMat[j] &= 0xf0;
//#endif
//        }
//        memcpy(matsRef[i]->vect, randMat, (*matsSecret)->dim_row * sizeof(BASE));
//
//        matTem = matSumOfRest;
//        matSumOfRest = add(matSumOfRest, matsRef[i]);
//        deMat(matTem);
//
//        /*  Refreshing  */
//        matTem = matsSecret[i];
//        matsSecret[i] = add(matsRef[i], matTem);
//        deMat(matTem);
//    }
//
//    free(randMat);
//
//#if DIM_A
//    matsRef[0] = multiply(matSumOfRest, matInvA);
//    deMat(matSumOfRest);
//
//      /* Refresh the first one, x_1 */
//    matTem = matsSecret[0];
//    matsSecret[0] = add(matsRef[0], matTem);
//    deMat(matTem);
//
//#else /* DIM_A == 0 */
//      /* Refresh the first one, x_1 */
//    matTem = matsSecret[0];
//    matsSecret[0] = add(matSumOfRest, matTem);
//    deMat(matTem);
//    deMat(matSumOfRest);
//#endif
//
//
//    return RES_OK;
//}




/* Add operation, as same as XOR */
Mat **addWithMask(
        const Mat **matEX,
        const Mat **matEY
        )
{
    if (matEX == NULL || matEY == NULL) return NULL;

    if (matEX[0]->dim_row != matEY[0]->dim_row ||
            matEX[0]->dim_col != matEY[0]->dim_col) return NULL;
    if (matEX[0]->vect == NULL || matEY[0]->vect == NULL) return NULL;

    /* Memory allocate */
    Mat **retMat = malloc(MASKD * sizeof(Mat *));
    if (retMat == NULL) return NULL;

    /* Calculate */
    int i;
    for (i = 0; i < MASKD; ++i){
        retMat[i] = add(matEX[i], matEY[i]);
    }

    return retMat;
}



/* Masked bitAnd operation  */
/* After operation, intermediate matrices will all be deallocated */



#if !DIVIDE && DIM_L == 16 && DIM_A
Mat **bitAndWithMask(
        const Mat **matEX,
        const Mat **matEY
        )
{
    if (matEX == NULL || matEY == NULL) return NULL;
    if (matEX[0]->dim_row != 1) return NULL;

    Mat *matTij[MASKD_SQURE], *matR[MASKD_SQURE];
    int i, j, k;
    int index;
    Mat *matTem;
    int slicesIndex;

    Mat **matEZ = (Mat **)malloc(MASKD * sizeof(Mat *));
    Mat *matEXPart[MASKD], *matEYPart[MASKD], *matEZPart[MASKD];

    for (i = 0; i != MASKD; ++i){
        matEZ[i] = newMat(1, DIM_L, NULL, 0x00);
        matEZPart[i] = newMat(1, DIM_L / DIVIDE_PARTS, NULL, 0x00);
        matEXPart[i] = newMat(1, DIM_L / DIVIDE_PARTS, NULL, 0x00);
        matEYPart[i] = newMat(1, DIM_L / DIVIDE_PARTS, NULL, 0x00);
    }
    for (slicesIndex = 0; slicesIndex != DIVIDE_PARTS; ++slicesIndex){
#if DIM_A == 4
        BYTE oddFlag = (BYTE)slicesIndex & 0x01;
#endif
        for (i = 0; i != MASKD; ++i){
#if DIM_A == 8
            matEXPart[i]->vect[0] = matEX[i]->vect[slicesIndex];
            matEYPart[i]->vect[0] = matEY[i]->vect[slicesIndex];
#elif DIM_A == 4
            int offset = (slicesIndex < 2) ? 0 : 1;
            matEXPart[i]->vect[0] = matEX[i]->vect[offset];
            matEYPart[i]->vect[0] = matEY[i]->vect[offset];

            if (oddFlag){//odd
                matEXPart[i]->vect[0] &= 0x0f;
                matEXPart[i]->vect[0] <<= 4;

                matEYPart[i]->vect[0] &= 0x0f;
                matEYPart[i]->vect[0] <<= 4;
            }
            else{
                matEXPart[i]->vect[0] &= 0xf0;
                matEYPart[i]->vect[0] &= 0xf0;
            }
#endif
        }
        /*  get matrix T  */
        for (i = 0; i != MASKD; ++i){
            for (j = 0; j != MASKD; ++j){
                index = i * MASKD + j;
#if DIM_A
                if (i == 0 && j == 0){
                    matTem = tensorProduct(matEXPart[i], matEYPart[j]);
                    matTij[index] = multiply(matTem, matHat);
                    deMat(matTem);
                }
                else if (i == 0){
                    matTem = tensorProduct(matEXPart[i], matEYPart[j]);
                    matTij[index] = multiply(matTem, matGrave);
                    deMat(matTem);
                }
                else if (j == 0){
                    matTem = tensorProduct(matEXPart[i], matEYPart[j]);
                    matTij[index] = multiply(matTem, matAcute);
                    deMat(matTem);
                }
                else {/*  i != 0 && j != 0  */
                    matTij[index] = bitAnd(matEXPart[i], matEYPart[j]);
                }

#else  /* DIM_A == 0 */
                matTij[index] = bitAnd(matEXPart[i], matEYPart[j]);
#endif /* DIM_A */
            }
        }
        int btsOfRow = bytesOfRow(matEXPart[0]->dim_col);
        int btsOfMat = matEXPart[0]->dim_row * btsOfRow;
        BASE *randMat = (BASE*)malloc(btsOfMat * sizeof(BASE));
        if (randMat == NULL) return NULL;

        /*  get matrix R  */
        for (i = 0; i < MASKD; ++i){

            /* get R(i,i) */
            index = i * MASKD + i;
            matR[index] = newMat(matEXPart[0]->dim_row, matEXPart[0]->dim_col, NULL, 0x00);
            memcpy(matR[index]->vect, matTij[index]->vect, btsOfMat * sizeof(BASE));
            deMat(matTij[index]);

            for (j = i + 1; j < MASKD; ++j){

                /* get R(i,j) through random generating */
                index = i * MASKD + j;
                matR[index] = newMat(matEXPart[0]->dim_row, matEXPart[0]->dim_col, NULL, 0x00);
                if (matR[index] == NULL) return NULL;

                for (k = 0; k < btsOfMat; ++k){
                    randMat[k] = (BASE)rand();

                    if (matEXPart[0]->dim_col == 4) randMat[k] &= 0xf0;

                }

                memcpy(matR[index]->vect, randMat, btsOfMat * sizeof(BASE));

                /* get R(j,i)  */
                matTem = add(matR[index], matTij[index]);
                deMat(matTij[index]);

                index = j * MASKD + i;
                matR[index] = add(matTij[index], matTem);

                deMat(matTij[index]);
                deMat(matTem);
#if DIM_A
                /*  get R(0,j) */
                if (i == 0){
                    index = j;
                    matTem = matR[index];
                    matR[index] = multiply(matTem, matA);
                    deMat(matTem);
                }
#endif  /* DIM_A != 0 */
            }
        }

        for (i = 0; i != MASKD; ++i){
            matEZPart[i] = matR[i];
            for (j = 1; j != MASKD; ++j){
                index = j * MASKD + i;
                matTem = matEZPart[i];
                matEZPart[i] = add(matTem, matR[index]);

                deMat(matTem);
                deMat(matR[index]);
            }
        }
        __noop;
        for (i = 0; i != MASKD; ++i){
#if DIM_A == 8
            matEZ[i]->vect[slicesIndex] = matEZPart[i]->vect[0];
#elif DIM_A == 4
            int offset = (slicesIndex < 2) ? 0 : 1;
            if (oddFlag){//odd
                BYTE tem = matEZPart[i]->vect[0] >> 4;
                matEZ[i]->vect[offset] ^= (tem & 0x0f);
            }else matEZ[i]->vect[offset] ^= matEZPart[i]->vect[0];

#endif
        }

        free(randMat);

    }
    for (i = 0; i != MASKD; ++i){
        deMat(matEZPart[i]);
        deMat(matEXPart[i]);
        deMat(matEYPart[i]);
    }

    return matEZ;
}

#else
Mat **bitAndWithMask(
        const Mat **matEX,
        const Mat **matEY
        )
{
    if (matEX == NULL || matEY == NULL) return NULL;

    Mat **matEZ = (Mat **)malloc(MASKD * sizeof(Mat *));
    if (matEZ == NULL) return NULL;

    Mat **matTij = (Mat **)malloc(MASKD * MASKD * sizeof(Mat *));
    if (matTij == NULL) return NULL;

    Mat **matR = (Mat **)malloc(MASKD * MASKD * sizeof(Mat *));
    if (matR == NULL) return NULL;

    int i, j, k;
    int index;
    Mat *matTem;
    /*  get matrix T  */
    for (i = 0; i < MASKD; ++i){
        for (j = 0; j < MASKD; ++j){
            index = i * MASKD + j;
#if DIM_A
            if (i == 0 && j == 0){
                matTem = tensorProduct(matEX[i], matEY[j]);
                matTij[index] = multiply(matTem, matHat);
                deMat(matTem);
            }
            else if (i == 0){
                matTem = tensorProduct(matEX[i], matEY[j]);
                matTij[index] = multiply(matTem, matGrave);
                deMat(matTem);
            }
            else if (j == 0){
                matTem = tensorProduct(matEX[i], matEY[j]);
                matTij[index] = multiply(matTem, matAcute);
                deMat(matTem);
            }
            else {/*  i != 0 && j != 0  */
                matTij[index] = bitAnd(matEX[i], matEY[j]);
            }

#else  /* DIM_A == 0 */
            matTij[index] = bitAnd(matEX[i], matEY[j]);
#endif /* DIM_A */
        }
    }
    int btsOfRow = bytesOfRow(matEX[0]->dim_col);
    int btsOfMat = matEX[0]->dim_row * btsOfRow;
    BASE *randMat = (BASE*)malloc(btsOfMat * sizeof(BASE));
    if (randMat == NULL) return NULL;

    /*  get matrix R  */
    for (i = 0; i < MASKD; ++i){

        /* get R(i,i) */
        index = i * MASKD + i;
        matR[index] = newMat(matEX[0]->dim_row, matEX[0]->dim_col, NULL, 0x00);
        memcpy(matR[index]->vect, matTij[index]->vect, btsOfMat * sizeof(BASE));
        deMat(matTij[index]);

        for (j = i + 1; j < MASKD; ++j){

            /* get R(i,j) through random generating */
            index = i * MASKD + j;
            matR[index] = newMat(matEX[0]->dim_row, matEX[0]->dim_col, NULL, 0x00);
            if (matR[index] == NULL) return NULL;

            for (k = 0; k < btsOfMat; ++k){
                randMat[k] = (BASE)rand();

                if (matEX[0]->dim_col == 4) randMat[k] &= 0xf0;

            }

            memcpy(matR[index]->vect, randMat, btsOfMat * sizeof(BASE));

            /* get R(j,i)  */
            matTem = add(matR[index], matTij[index]);
            deMat(matTij[index]);

            index = j * MASKD + i;
            matR[index] = add(matTij[index], matTem);

            deMat(matTij[index]);
            deMat(matTem);
#if DIM_A
            /*  get R(0,j) */
            if (i == 0){
                index = j;
                matTem = matR[index];
                matR[index] = multiply(matTem, matA);
                deMat(matTem);
            }
#endif  /* DIM_A != 0 */
        }
    }

    for (i = 0; i < MASKD; ++i){
        matEZ[i] = matR[i];
        for (j = 1; j < MASKD; ++j){
            index = j * MASKD + i;
            matTem = matEZ[i];
            matEZ[i] = add(matTem, matR[index]);

            deMat(matTem);
            deMat(matR[index]);
        }
    }

    /* Deallocating*/
    free(randMat);
    free(matTij);
    free(matR);

    return matEZ;
}
#endif

///* Matrix multiplication (Original) */
//Mat *multiplyOriginal(
//              const Mat *matX,
//              const Mat *matY
//              )
//{
//    if (matX == NULL || matY == NULL) return NULL;
//    if (matX->dim_col != matY->dim_row) return NULL;
//    if (matX->vect == NULL || matY->vect == NULL) return NULL;
//
//    /* Memory allocation */
//    Mat *matRet = newMat(matX->dim_row, matY->dim_col, NULL, 0x00);
//    if (matRet == NULL) return NULL;
//
//    /* Preparing a transposed matrix first */
//    Mat *matYT = transpose(matY);
//
//    /* Multipling */
//    int row, col;
//    BASE *ptrOfMatRet, *ptrOfMatX, *ptrOfMatYT;
//    for (row = 0; row < matX->dim_row; ++row){
//        for (col = 0; col < matY->dim_col; ++col)
//        {
//            int bytesOfRX = bytesOfRow(matX->dim_col);
//            int bytesOfRR = bytesOfRow(matY->dim_col);
//
//            int cntsVect = col / LENGTH;
//            int offset = col % LENGTH;
//            ptrOfMatRet = matRet->vect + row * bytesOfRR + cntsVect;
//
//            int i;
//            BASE vectTem;
//            for (i = 0; i < bytesOfRX; ++i){
//                ptrOfMatX = matX->vect + row * bytesOfRX + i;
//                ptrOfMatYT = matYT->vect + col * bytesOfRX + i;
//                vectTem = (*ptrOfMatX) & (*ptrOfMatYT);
//                vectTem = sumOfByte(vectTem);
//                vectTem = vectTem >> offset;
//                (*ptrOfMatRet) ^= vectTem;
//            }
//        }
//    }
//
//    deMat(matYT);
//    return matRet;
//}

#endif /* MASK */


/*  Matrix Multiply(Transposed), (A, B) == A x B^T   */
Mat *multiply(
        const Mat *matX,
        const Mat *matY
        )
{
    if (matX == NULL || matY == NULL) return NULL;
    if (matX->dim_col != matY->dim_col) return NULL;
    if (matX->vect == NULL || matY->vect == NULL) return NULL;

    /* Memory allocation */
    Mat *matRet = newMat(matX->dim_row, matY->dim_row, NULL, 0x00);
    if (matRet == NULL) return NULL;


    /* Multipling */
    int row, col;
    BASE *ptrOfMatRet, *ptrOfMatX, *ptrOfMatY;
    for (row = 0; row < matX->dim_row; ++row){
        for (col = 0; col < matY->dim_row; ++col)
        {
            int bytesOfRX = bytesOfRow(matX->dim_col);
            int bytesOfRR = bytesOfRow(matY->dim_row);

            int cntsVect = col / LENGTH;
            int offset = col % LENGTH;
            ptrOfMatRet = matRet->vect + row * bytesOfRR + cntsVect;

            int i;
            BASE vectTem;
            for (i = 0; i < bytesOfRX; ++i){
                ptrOfMatX = matX->vect + row * bytesOfRX + i;
                ptrOfMatY = matY->vect + col * bytesOfRX + i;

                // vectTem = sumOfByte( (*ptrOfMatX) & (*ptrOfMatY) );
				// vectTem >>= offset;
				vectTem = (BASE)_mm_popcnt_u32((*ptrOfMatX) & (*ptrOfMatY)) & 0x01;
				vectTem <<= (LENGTH -1 - offset);

                (*ptrOfMatRet) ^= vectTem;
            }
        }
    }

    return matRet;
}


/* Simple Add operation, as same as XOR */
Mat *add(
        const Mat *matX,
        const Mat *matY
        )
{
    if (matX == NULL || matY == NULL) return NULL;

    if (matX->dim_row != matY->dim_row ||
            matX->dim_col != matY->dim_col) return NULL;
    if (matX->vect == NULL || matY->vect == NULL) return NULL;

    /* Memory allocate */
    int col = matX->dim_col;
    int row = matX->dim_row;

    Mat *retMat = newMat(row, col, NULL, 0x00);

    if (retMat == NULL) return NULL;

    /* Calculate */
    int i;
    int bytesOfMat = row * bytesOfRow(col);
    for (i = 0; i != bytesOfMat; ++i){
        retMat->vect[i] = matX->vect[i] ^ matY->vect[i];
    }

    return retMat;
}



/* Simple bitAnd operation , as same as AND */
Mat *bitAnd(
        const Mat *matX,
        const Mat *matY
        )
{
    if (matX == NULL || matY == NULL) return NULL;

    if (matX->dim_row != matY->dim_row ||
            matX->dim_col != matY->dim_col) return NULL;
    if (matX->vect == NULL || matY->vect == NULL) return NULL;

    /* memory allocate */
    Mat *retMat = newMat(matX->dim_row, matX->dim_col, NULL, 0x00);
    if (retMat == NULL) return NULL;

    /* calculate */
    int i;
    int bytesOfMat = matX->dim_row * bytesOfRow(matX->dim_col);
    for (i = 0; i != bytesOfMat; ++i){
        retMat->vect[i] = matX->vect[i] & matY->vect[i];
    }

    return retMat;
}




/* Split a matrix to n parts(sub-mats) through the dimension r */
Mat **split(
        const Mat *matOrig,
        int n,
        int r
        /* r can only be '1' or '2'
         * '1' : split  through row-dimension
         * '2' :     ...        col-dimension,
         */
        )
{
    if (matOrig == NULL) return NULL;
    int index, subR, subC, bytesOfSubMat;

    Mat **retMats = (Mat **)malloc(n * sizeof(Mat *));
    if (retMats == NULL) return NULL;

    BASE *ptrOfSubMat;
    if (r == 1){
        if (matOrig->dim_row % n) return NULL;
        subR = matOrig->dim_row / n;
        subC = matOrig->dim_col;
        bytesOfSubMat = subR * bytesOfRow(matOrig->dim_col);
        ptrOfSubMat = matOrig->vect;
        for (index = 0; index < n; ++index){
            retMats[index] = newMat(subR, subC, NULL, 0x00);
            if (retMats[index] == NULL) return NULL;
            memcpy(retMats[index]->vect, ptrOfSubMat, bytesOfSubMat);
            ptrOfSubMat += bytesOfSubMat;
        }
        ptrOfSubMat = NULL;
    }
    else if (r == 2){
        /* other cases are To Be Implemented */

        subR = matOrig->dim_row;
        subC = matOrig->dim_col / n;

        bytesOfSubMat = subR * bytesOfRow(subC);

        Mat *transOrig = transpose(matOrig);
        Mat **transRes = split(transOrig, n, 1);
        deMat(transOrig);
        for (index = 0; index < n; ++index){
            retMats[index] = transpose(transRes[index]);
            deMat(transRes[index]);
        }
    }
    else return NULL;
    return retMats;
}



/* Catenate n mats through dimension r */
Mat *cat(
        const Mat **mats,
        int n,
        int r
        )
{
    if (mats == NULL || *mats == NULL) return NULL;

    int index, subR, subC, bytesOfSubMat;
    BASE *ptrOfBigMat;
    Mat *retMat;
    subR = (*mats)->dim_row;
    subC = (*mats)->dim_col;

    if (r == 1){

        retMat = newMat(subR * n, subC, NULL, 0x00);
        if (retMat == NULL) return NULL;
        bytesOfSubMat = subR * bytesOfRow(subC);
        ptrOfBigMat = retMat->vect;
        for (index = 0; index < n; ++index){
            if (mats[index] == NULL) return NULL;
            memcpy(ptrOfBigMat, mats[index]->vect, bytesOfSubMat);
            ptrOfBigMat += bytesOfSubMat;
        }
        ptrOfBigMat = NULL;
    }
    else if (r == 2){
        /* other cases are To Be Implemented */
        Mat **transMats = (Mat **)malloc(n * sizeof(Mat *));
        if (transMats == NULL) return NULL;
        for (index = 0; index < n; ++index){
            transMats[index] = transpose(mats[index]);
        }

        Mat *transRes = cat(transMats, n, 1);
        retMat = transpose(transRes);
        deMat(transRes);
        deMats(transMats, n);
    }
    else return NULL;

    return retMat;

}


/* =========================================== */
/*                The   End                    */
/* =========================================== */
