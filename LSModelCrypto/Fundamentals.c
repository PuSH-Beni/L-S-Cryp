//
//  Fundamentals.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//


#include "Fundamentals.h"

/*=====================*/
/* private functions */
/*=====================*/


/* caculate the num of bytes in one row */
static
int bytesOfRow(
    int col
)
{
    // bytes of each row :: if dim_col < 8, allocate a byte as well
    int countOfBytes = col >> 3; // col / 8
    if(countOfBytes == 0)countOfBytes = 1;
    return countOfBytes;
}


/* naive 'and' operation, return a byte with posbase as res */
static
BYTE and(
    const BYTE byteBase,
    int posBase,
    const BYTE byteAdjust,
    int posAdjust
)
{
    BYTE res;
    if(posBase == posAdjust) res = byteBase & byteAdjust;
    else if(posBase > posAdjust){
        res = byteAdjust;
        res = res >> (posBase - posAdjust);
        res &= byteBase;
    }else{
        res = byteAdjust;
        res = res << (posAdjust - posBase);
        res &= byteBase;
    }
    BYTE tem = 0x80;// 1000_0000 for binary
    res &= (tem >> posBase);
    return res;
}




/*=====================*/
/* pubulic functions */
/*=====================*/

/* construce a new Mat instance */
Mat *newMat(
    int dim_row,
    int dim_col,
    BYTE flags /* flags:
                *   '0': norm
                *   '1': unit matrix or indentit matrix
                *   '2': error
                */
)
{
    if(dim_row <= 0 || dim_col <= 0)return NULL;
    Mat *retMat = (Mat *)malloc(sizeof(Mat));
    if(retMat == NULL) return NULL;

    int countOfBytes = bytesOfRow(dim_col) * dim_row;
    retMat->vect = (BYTE *)malloc(countOfBytes * sizeof(BYTE));
    memset(retMat->vect,0,countOfBytes * sizeof(BYTE));
    if(retMat->vect == NULL)
    {
        free(retMat);
        return NULL;
    }
    retMat->dim_row = dim_row;
    retMat->dim_col = dim_col;
	retMat->flags = flags;
    return retMat;
}

/* deconstruct a Mat instance */
void deMat(
    Mat *matrix
)
{
    if(matrix != NULL)
    {
        free(matrix->vect);
        free(matrix);
    }
}

/* add operation, as same as XOR*/
Mat *add(
    const Mat *matA,
    const Mat *matB
)
{
    if(matA == NULL || matB == NULL) return NULL;

    if(matA->dim_row != matB->dim_row ||
            matA->dim_col != matB->dim_col) return NULL;
    if(matA->vect == NULL || matB->vect == NULL) return NULL;

    Mat *retMat = newMat(matA->dim_row, matA->dim_col, 0);
    if(retMat == NULL) return NULL;
    int countOfBytes = bytesOfRow(matB->dim_col) * matB->dim_row;

    if(matB->flags == 1){
        memmove(retMat->vect, matA->vect, countOfBytes);
        retMat->flags = matA->flags;
    }
    else if(matA->flags == 1){
        memmove(retMat->vect, matB->vect, countOfBytes);
        retMat->flags = matB->flags;
    }
    else{
        *(retMat->vect) = *(matA->vect) ^ *(matB->vect);
    }
    return retMat;
}


/* bitAnd operation without mask, is just naive 'and' logic */
Mat *bitAnd(
    const Mat *matA,
    const Mat *matB
)
{
    if(matA == NULL || matB == NULL) return NULL;
    if(matA->dim_row != matB->dim_row ||
            matA->dim_col != matB->dim_col) return NULL;
    if(matA->vect == NULL || matB->vect == NULL) return NULL;

    int countOfBytes = bytesOfRow(matB->dim_col) * matB->dim_row;

    Mat *retMat = newMat(matA->dim_row, matA->dim_col, 0);
    if(retMat == NULL) return NULL;

    if(matB->flags == 1){
        memmove(retMat->vect, matA->vect, countOfBytes);
        retMat->flags = matA->flags;
    }
    else if(matA->flags == 1){
        memmove(retMat->vect, matB->vect, countOfBytes);
        retMat->flags = matB->flags;
    }
    else{
        *(retMat->vect) = *(matA->vect) & *(matB->vect);
    }
	return retMat;
}


/* matrix multiply */
Mat *multiply(
    const Mat *matA,
    const Mat *matB
)
{
    if(matA == NULL || matB == NULL) return NULL;
    if(matA->dim_col != matB->dim_row) return NULL;
    if(matA->vect == NULL || matB->vect == NULL) return NULL;

    int countOfBytes;
    Mat *retMat = newMat(matA->dim_row, matB->dim_col, 0);
    if(retMat == NULL) return NULL;

    int row, col, it;

    /* define the step, which is the flag to forward the byte */
    int stepMatB = matB->dim_col < 8 ? matB->dim_col : 8;
    int stepMatA = matA->dim_col < 8 ? matA->dim_col : 8;
    int bytesOfRowA = bytesOfRow(matA->dim_col);
    int bytesOfRowB = bytesOfRow(matB->dim_col);
    int posMatA, posMatB;
    BYTE *ptrOfRet = retMat->vect;
    BYTE *ptrOfMatA = matA->vect, *ptrOfMatB = matB->vect;

    /* element from matrix B and matrix Ret */
    BYTE elemR;
    for(row = 0; row < matA->dim_row; ++row){
        ptrOfMatA = matA->vect + row * bytesOfRowA;
        for(col = 0; col < matB->dim_col; ++col){
            if(col >= stepMatB){
                ptrOfRet += 1;
                ptrOfMatB += 1;
                posMatB = col - stepMatB;
            }else{
                posMatB = col;
            }
            elemR = 0;
            for(it = 0; it < matA->dim_col; ++it){
                ptrOfMatB = matB->vect + it * bytesOfRowB;
                if(col >= stepMatB){
                    ptrOfMatB  += 1;
                }
                if(it >= stepMatA){
                    ptrOfMatA += 1;
                    posMatA = it - stepMatA;
                }else{
                    posMatA = it;
                }
                elemR ^= and(*ptrOfMatB, posMatB, *ptrOfMatA, posMatA);
            }
            *ptrOfRet ^= elemR;
        }
        if(row < (matA->dim_row - 1))ptrOfRet += 1;
        else ptrOfRet = NULL;
    }
    return retMat;
}

/* split a matrix to n parts(sub-mat) through the dimension r */
Mat **split(
    const Mat *matO,
    const int n,
    const int r /* r can only be '1' or '2'
           * '1' :: split  through row-dimension
           * '2' ::     ...        col-dimension
           */
)
{
    if(matO == NULL) return NULL;
    int index, subR, subC, bytesOfSubMat;

    Mat **retMats = (Mat **)malloc(n*sizeof(Mat *));
    if(retMats == NULL) return NULL;

    BYTE *ptrOfSubMat;
    if(r == 1){
        if(matO->dim_row % n) return NULL;
        subR = matO->dim_row / n;
        subC = matO->dim_col;
        bytesOfSubMat = subR * bytesOfRow(matO->dim_col);
        ptrOfSubMat = matO->vect;
        for(index = 0; index < n; ++index){
            retMats[index] = newMat(subR, subC, 0);
            if(retMats[index] == NULL) return NULL;
            memmove(retMats[index]->vect, ptrOfSubMat, bytesOfSubMat);
            ptrOfSubMat += bytesOfSubMat;
        }
        ptrOfSubMat = NULL;
    }else if(r == 2){

        /* To Be Implemented */
        return NULL;

    }else return NULL;
    return retMats;
}


/* catenate n mats through dimension r */
Mat *cat(
    const Mat **mats,
    const int n,
    const int r
)
{
    if(mats == NULL || *mats == NULL ) return NULL;

    int index, subR, subC, bytesOfSubMat;
    BYTE *ptrOfBigMat;
    Mat *retMat;
    subR = (*mats)->dim_row;
    subC = (*mats)->dim_col;

    if(r == 1){

        retMat = newMat(subR * n, subC, 0);
        if(retMat == NULL) return NULL;
        bytesOfSubMat = subR * bytesOfRow(subC);
        ptrOfBigMat = retMat->vect;
        for(index = 0; index < n; ++index){
            if(mats[index] == NULL) return RES_INVALID_POINTER;
            memmove(ptrOfBigMat, mats[index]->vect, bytesOfSubMat);
            ptrOfBigMat += bytesOfSubMat;
        }
        ptrOfBigMat = NULL;
    }else if(r == 2){

        /* To Be Implemented */
        return NULL;

    }else return NULL;

    return retMat;

}
