//
//  Fundamentals.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//


#include "Fundamentals.h"

/*=========================================================*/
/*   Private   Functions      */
/*=========================================================*/


/* Caculate the num of bytes in one row */
static
int bytesOfRow(
int col
)
{
	// bytes of each row :: if dim_col < 8, allocate a byte as well
	int countOfBytes = col >> 3; // col / 8
	if (countOfBytes == 0)countOfBytes = 1;
	return countOfBytes;
}


/* Another AND operation for two bytes, return a byte with posbase as res */
static
BYTE and(
const BYTE byteBase,
int posBase,
const BYTE byteAdjust,
int posAdjust
)
{
	BYTE res;
	if (posBase == posAdjust) res = byteBase & byteAdjust;
	else if (posBase > posAdjust){
		res = byteAdjust;
		res = res >> (posBase - posAdjust);
		res &= byteBase;
	}
	else{
		res = byteAdjust;
		res = res << (posAdjust - posBase);
		res &= byteBase;
	}
	BYTE tem = 0x80;// 1000_0000 for binary
	res &= (tem >> posBase);
	return res;
}



/* Transpose a vector from (j,i) to (i,j) */
VECT transVector(
	VECT orig,
	int i,
	int j
	)
{
	VECT unit = IDENT >> i;
	VECT tem = unit & orig;

	if (i > j) tem = tem << (i - j);
	else if (i < j) tem = tem >> (j - i);
	else return tem;

	return tem;
}


/* Get sum of a row */
VECT sumFromRow(
	VECT row,
	int cols
	)
{
	int i;
	VECT ret = ZEROV;
	for (i = 0; i < cols; ++i){
		ret = ret ^ (row << i);
	}
	ret &= IDENT;
	return ret;
}


/*=========================================================*/
/*   Public   Functions      */
/*=========================================================*/

/* Construct a new Mat instance */
Mat *newMat(
	int dim_row,
	int dim_col,
	VECT *addr,
	BYTE flags /*   flags :
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
	if (addr == NULL){
		retMat->vect = (VECT *)malloc(dim_row * sizeof(VECT));
		memset(retMat->vect, 0, dim_row * sizeof(VECT));
		retMat->flags = flags;
	}
	else{
		retMat->vect = addr;
		retMat->flags = 0x03;
	}
	
	
	if (retMat->vect == NULL)
	{
		free(retMat);
		return NULL;
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


/* Transpositon */
Mat *transpose(
	const Mat *matA
	)
{
	if (matA == NULL) return NULL;
	int row, col;
	VECT rowToTrans, rowBeTransed;

	/* Memory allocated */
	Mat *matRet = newMat(matA->dim_col, matA->dim_row, NULL, 0x00);
	if (matRet == NULL) return NULL;

	/* Transposing */
	for (row = 0; row < matA->dim_col; ++row){
		matRet->vect[row] = ZEROV;
		for (col = 0; col < matA->dim_row; ++col){
			rowToTrans = matA->vect[col];
			rowBeTransed = transVector(rowToTrans, row, col);
			matRet->vect[row] ^= rowBeTransed;
		}
	}
	return matRet;

}



/* Add operation, as same as XOR */
Mat *add(
	const Mat *matA,
	const Mat *matB
	)
{
	if (matA == NULL || matB == NULL) return NULL;

	if (matA->dim_row != matB->dim_row ||
		matA->dim_col != matB->dim_col) return NULL;
	if (matA->vect == NULL || matB->vect == NULL) return NULL;

	/* memory allocate */
	Mat *retMat = newMat(matA->dim_row, matA->dim_col, NULL, 0x00);
	if (retMat == NULL) return NULL;

	/* calculate */
	int i;
	for (i = 0; i < matA->dim_row; ++i){
		retMat->vect[i] = matA->vect[i] ^ matB->vect[i];
	}

	return retMat;
}


/* BitAnd operation , as same as AND */
Mat *bitAnd(
	const Mat *matA,
	const Mat *matB
	)
{
	if (matA == NULL || matB == NULL) return NULL;

	if (matA->dim_row != matB->dim_row ||
		matA->dim_col != matB->dim_col) return NULL;
	if (matA->vect == NULL || matB->vect == NULL) return NULL;

	/* memory allocate */
	Mat *retMat = newMat(matA->dim_row, matA->dim_col, NULL, 0x00);
	if (retMat == NULL) return NULL;

	/* calculate */
	int i;
	for (i = 0; i < matA->dim_row; ++i){
		retMat->vect[i] = matA->vect[i] & matB->vect[i];
	}

	return retMat;
}


/* Matrix multiplication */
Mat *multiply(
	const Mat *matA,
	const Mat *matB
	)
{
	if (matA == NULL || matB == NULL) return NULL;
	if (matA->dim_col != matB->dim_row) return NULL;
	if (matA->vect == NULL || matB->vect == NULL) return NULL;

	/* Memory allocation */
	Mat *matRet = newMat(matA->dim_row, matB->dim_col, NULL, 0x00);
	if (matRet == NULL) return NULL;

	/* Preparing a transposed matrix first */
	Mat *matBT = transpose(matB);

	/* Multipling */
	int row, col;
	for (row = 0; row < matA->dim_row; ++row){
		VECT rowBeCaled = ZEROV;
		for (col = 0; col < matB->dim_col; ++col)
		{
			VECT rowTem;
			rowTem = matA->vect[row] & matBT->vect[col];
			rowTem = sumFromRow(rowTem, matA->dim_col);
			rowTem = rowTem >> col;
			rowBeCaled ^= rowTem;
		}
		matRet->vect[row] = rowBeCaled;
	}

	return matRet;
}


	

/* Split a matrix to n parts(sub-mats) through the dimension r */
Mat **split(
	const Mat *matO,
	const int n,
	const int r 
		  /* r can only be '1' or '2'
		   * '1' :: split  through row-dimension
		   * '2' ::     ...        col-dimension
		   */
		   )
{
	if (matO == NULL) return NULL;
	int index, subR, subC, bytesOfSubMat;

	Mat **retMats = (Mat **)malloc(n * sizeof(Mat *));
	if (retMats == NULL) return NULL;

	VECT *ptrOfSubMat;
	if (r == 1){
		if (matO->dim_row % n) return NULL;
		subR = matO->dim_row / n;
		subC = matO->dim_col;
		bytesOfSubMat = subR * bytesOfRow(matO->dim_col);
		ptrOfSubMat = matO->vect;
		for (index = 0; index < n; ++index){
			retMats[index] = newMat(subR, subC, NULL, 0x00);
			if (retMats[index] == NULL) return NULL;
			memmove(retMats[index]->vect, ptrOfSubMat, bytesOfSubMat);
			ptrOfSubMat += bytesOfSubMat;
		}
		ptrOfSubMat = NULL;
	}
	else if (r == 2){

		/* To Be Implemented */
		return NULL;

	}
	else return NULL;


	return retMats;
}


/* Catenate n mats through dimension r */
Mat *cat(
	const Mat **mats,
	const int n,
	const int r
	)
{
	if (mats == NULL || *mats == NULL) return NULL;

	int index, subR, subC, bytesOfSubMat;
	VECT *ptrOfBigMat;
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
			memmove(ptrOfBigMat, mats[index]->vect, bytesOfSubMat);
			ptrOfBigMat += bytesOfSubMat;
		}
		ptrOfBigMat = NULL;
	}
	else if (r == 2){

		/* To Be Implemented */
		return NULL;

	}
	else return NULL;

	return retMat;

}
