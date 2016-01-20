//
//  Fundamentals.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//


#include "Fundamentals.h"



/*=========================================================*/
/*   Declarition     */
/*=========================================================*/
#if !TEST
static
#endif
Mat *matA, *matInvA, *matTransA, *matHat, *matGrave, *matAcute;



/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/


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
#if LENG8
	bytes = col >> 3; // col / LENGTH
#elif LENG4
	bytes= col >> 2; // col / LENGTH
#endif
	bytes = (col % LENGTH == 0) ? bytes : bytes + 1;
	return bytes;
}


/* Another AND operation for two bytes, return a byte with posbase as res */
/*static
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
*/


/* Shift bit form j -> i */
#if !TEST
static
#endif 
BYTE shiftBit(
const BYTE orig,
const int i,
const int j
)
{
	BYTE unit = IDENT >> j;
	BYTE tem = unit & orig;

	if (j > i) tem = tem << (j - i);
	else if (j < i) tem = tem >> (i - j);
	else return tem;

	return tem;
}



/* Get sum of a Vect */
#if !TEST
static
#endif 
BYTE sumFromVect(
BYTE Vect
)
{

	BYTE ret = 0;
	while (Vect > 0){
		ret ^= IDENT;
		Vect &= (Vect - 1);
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
	//srand((WORD)time(NULL));

	int *nums = (int *)malloc(LENGTH * sizeof(int));
	if (nums == NULL) return NULL;
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
const Mat *matO
)
{
	if (matO == NULL) return NULL;

	/* Memory allocated */
	Mat *matRet = newMat(matO->dim_col, matO->dim_row, NULL, 0x00);
	if (matRet == NULL) return NULL;

	int bytesOfRO = bytesOfRow(matO->dim_col);
	int bytesOfRR = bytesOfRow(matO->dim_row);

	/* Transposing */
	int colO, rowO;
	int cntbytesO, cntbytesR;
	int offsetO, offsetR;
	BYTE vectTransed;
	BYTE *ptrOfVectRet, *ptrOfVectO;
	for (colO = 0; colO < matO->dim_col; ++colO){
		cntbytesO = colO / LENGTH;
		offsetO = colO % LENGTH;
		for (rowO = 0; rowO < matO->dim_row; ++rowO){
			/* the bit (rowO, colO) */
			cntbytesR = rowO / LENGTH;
			offsetR = rowO % LENGTH;

			ptrOfVectO = matO->vect + bytesOfRO * rowO + cntbytesO;
			ptrOfVectRet = matRet->vect + bytesOfRR * colO + cntbytesR;

			vectTransed = shiftBit(*ptrOfVectO, offsetR, offsetO);
			(*ptrOfVectRet) ^= vectTransed;
		}
	}
	return matRet;

}





/* ============================================================== */
/*      Functions For Masked Model                  */
/* ============================================================== */
#if MASK
/* Generate LENGTH x LENGTH random matrices A, inv(A), trans(A)  */

#if !TEST
static
#endif 
Mat **genRandMat(
)
{

	/* Memory allocated */
	Mat **retMats = (Mat **)malloc(3 * sizeof(Mat *));
	if (retMats == NULL) return NULL;
	BYTE vecId[] = INDENT_MAT;
	Mat *matE = newMat(LENGTH, LENGTH, vecId, 0x03);
	Mat *matATem = newMat(LENGTH, LENGTH, NULL, 0x00);
	Mat  *matTransP = NULL, *matTem;

	/*  get matrix A and inv(A) temporarily */
	int *randRow = randOrder();
	int i;
	for (i = 0; i < LENGTH; ++i) matATem->vect[i] = matE->vect[i];
	Mat *matInvATem = transpose(matATem);
	Mat *matTransATem = transpose(matATem);
	deMat(matATem);

	/* get matrix P */
	//srand((WORD)time(NULL));
	for (i = 0; i < LENGTH; ++i){
		BYTE rowToAdd = (BYTE)rand();
		BYTE zeroToSet = IDENT >> i;
		rowToAdd &= (~zeroToSet);

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
	free(randRow);
	deMat(matTransP);
	deMat(matE);

	return retMats;
}




/* Generate hatA, graveA and acuteA */
/* All of them are n x n^2 */
#if !TEST
static
#endif 
Mat *hatA(
)
{
	if (matTransA == NULL || matInvA == NULL) return NULL;

	Mat *matRight = newMat(LENGTH * LENGTH, LENGTH, NULL, 0x00);
	if (matRight == NULL) return NULL;

	/* get transposition of the right matrix: (E' x (A tp B))^T */
	int i, j;
	BYTE *ptrOfMatATi;
	BYTE *ptrOfMatATj, *ptrOfMatR;
	for (i = 0; i < LENGTH; ++i){
		ptrOfMatATi = matTransA->vect + i;
		for (j = 0; j < LENGTH; ++j){
			ptrOfMatATj = matTransA->vect + j;
			ptrOfMatR = matRight->vect + (i * LENGTH + j);
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
Mat *graveA(
)
{
	if (matTransA == NULL || matInvA == NULL) return NULL;

	Mat *matRight = newMat(LENGTH * LENGTH, LENGTH, NULL, 0x00);
	if (matRight == NULL) return NULL;


	int i, j;
	BYTE *ptrOfMatATj, *ptrOfMatR, *ptrOfMatI;
	BYTE vecId[] = INDENT_MAT;
	Mat *matIdent = newMat(LENGTH, LENGTH, vecId, 0x03);
	if (matIdent == NULL) return NULL;


	/* get transposition of the right matrix: (E' x (A tp E))^T */
	for (i = 0; i < LENGTH; ++i){
		ptrOfMatI = matIdent->vect + i;
		for (j = 0; j < LENGTH; ++j){
			ptrOfMatATj = matTransA->vect + j;
			ptrOfMatR = matRight->vect + (i * LENGTH + j);
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
Mat *acuteA(
)
{

	if (matTransA == NULL || matInvA == NULL) return NULL;

	Mat *matRight = newMat(LENGTH * LENGTH, LENGTH, NULL, 0x00);
	if (matRight == NULL) return NULL;


	int i, j;
	BYTE *ptrOfMatATi, *ptrOfMatR, *ptrOfMatI;
	/* Identity matrix */
	BYTE vecId[] = INDENT_MAT;
	Mat *matIdent = newMat(LENGTH, LENGTH, vecId, 0x03);
	if (matIdent == NULL) return NULL;


	/* get transposition of the right matrix: (E' x (E tp A))^T */
	for (i = 0; i < LENGTH; ++i){
		ptrOfMatATi = matTransA->vect + i;
		for (j = 0; j < LENGTH; ++j){
			ptrOfMatI = matIdent->vect + j;
			ptrOfMatR = matRight->vect + i * LENGTH + j;
			(*ptrOfMatR) = (*ptrOfMatI) & (*ptrOfMatATi);
		}
	}

	/* get \Acute{A} now */
	Mat *matRet = multiply(matInvA, matRight);

	deMat(matRight);
	deMat(matIdent);
	return matRet;
}


/* Tenser Product for two vectors */
/* Considering 'matX->dim_col' is '8' or '4', a row of the matrix will be no longer than a BYTE */
/* Return a n x n^2 matrix */
#if !TEST
static
#endif 
Mat *tenserProduct(
const Mat *matX,
const Mat *matY
)
{
	if (matX == NULL || matY == NULL) return NULL;
	if ((matX->dim_row != matY->dim_row) ||
		(matX->dim_col != matY->dim_col))return NULL;

	if (matX->dim_col != LENGTH) return NULL;

	int i, j;
	Mat *matR = newMat(matX->dim_row, LENGTH * LENGTH, NULL, 0x00);
	if (matR == NULL) return NULL;

	for (i = 0; i < matX->dim_row; ++i){
		BYTE ident = IDENT;
		for (j = 0; j < LENGTH; ++j){
			matR->vect[i * LENGTH + j] = (matY->vect[i] & ident) ? matX->vect[i] : 0x00;
			ident = ident >> 1;
		}
	}

	return matR;
}

#endif


/*=========================================================*/
/*   MARK: Public   Functions      */
/*=========================================================*/

/* Construct a new Mat instance */
Mat *newMat(
	int dim_row,
	int dim_col,
	BYTE *addr,
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
	int bytesOfMat = dim_row * bytesOfRow(dim_col);
	if (addr == NULL){
		retMat->vect = (BYTE *)malloc(bytesOfMat * sizeof(BYTE));
		if (retMat->vect == NULL)
		{
			free(retMat);
			return NULL;
		}
		memset(retMat->vect, 0, bytesOfMat * sizeof(BYTE));
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






/* ============================================================== */
/*      Functions For Masked Model                  */
/* ============================================================== */

#if MASK

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


	free(mats);
}


/* Encode(Mask) the plain text */
Mat **encode(
	const Mat *matPlain
	)
{
	if (matInvA == NULL || matPlain == NULL) return NULL;
	Mat **matsRet = (Mat **)malloc(MASKD * sizeof(Mat *));
	if (matsRet == NULL) return NULL;

	int i, j;
	BYTE *randMat = (BYTE *)malloc(matPlain->dim_row  * sizeof(BYTE));
	if (randMat == NULL) return NULL;

	BYTE zeroV[] = ZERO_MAT;
	Mat *matSumOfRest = newMat(matPlain->dim_row, matPlain->dim_col, zeroV, 0x03);
	if (matSumOfRest == NULL) return NULL;
	
	Mat *matTem = NULL;

	for (i = 1; i < MASKD; ++i){

		matsRet[i] = newMat(matPlain->dim_row, matPlain->dim_col, NULL, 0x00);
		if (matsRet[i] == NULL) return NULL;

		for (j = 0; j < matPlain->dim_row; ++j){
			randMat[j] = (BYTE)rand();
		}
		memmove(matsRet[i]->vect, randMat, matPlain->dim_row * sizeof(BYTE));

		if (matsRet[i] == NULL) return NULL;

		matTem = matSumOfRest;
		matSumOfRest = add(matSumOfRest, matsRet[i]);
		deMat(matTem);

	}

	free(randMat);
	matTem = add(matSumOfRest, matPlain);
	deMat(matSumOfRest);


	/*  = (x^T) x inv(A) ^ T */
	matsRet[0] = multiply(matTem, matInvA);
	deMat(matTem);

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
	Mat *matSum = multiply(matsSecret[0], matA);
	deMat(matsSecret[0]);
	for (i = 1; i < MASKD; ++i){
		matOrig = matSum;
		matSum = add(matSum, matsSecret[i]);
		deMat(matOrig);
		deMat(matsSecret[i]);
	}
	free(matsSecret);
	return matSum;
}

/* Refresh the old mask with a new mask */
Res refreshing(
	Mat **matsSecret
	)
{
	if (matInvA == NULL || *matsSecret == NULL) return RES_INVALID_POINTER;
	Mat **matsRef = (Mat **)malloc(MASKD * sizeof(Mat *));
	if (matsRef == NULL) return RES_INVALID_POINTER;

	int i, j;
	BYTE *randMat = (BYTE *)malloc((*matsSecret)->dim_row * sizeof(BYTE));
	if (randMat == NULL) return RES_INVALID_POINTER;

	/* a ZERO matrix */
	BYTE zeroV[] = ZERO_MAT;

	Mat *matSumOfRest = newMat((*matsSecret)->dim_row, (*matsSecret)->dim_col, zeroV, 0x03);
	if (matSumOfRest == NULL) return RES_INVALID_POINTER;

	Mat *matTem = NULL;

	for (i = 1; i < MASKD; ++i){
		matsRef[i] = newMat(LENGTH, (*matsSecret)->dim_col, NULL, 0x00);
		if (matsRef[i] == NULL) return RES_INVALID_POINTER;

		for (j = 0; j < (*matsSecret)->dim_row; ++j){
			randMat[j] = (BYTE)rand();
		}
		memmove(matsRef[i]->vect, randMat, (*matsSecret)->dim_row * sizeof(BYTE));

		matTem = matSumOfRest;
		matSumOfRest = add(matSumOfRest, matsRef[i]);
		deMat(matTem);

		/*  Refreshing  */
		matTem = matsSecret[i];
		matsSecret[i] = add(matsRef[i], matTem);
		deMat(matTem);
	}

	free(randMat);


	matsRef[0] = multiply(matSumOfRest, matInvA);
	deMat(matSumOfRest);

	/* Refresh the first one, x_1 */
	matTem = matsSecret[0];
	matsSecret[0] = add(matsRef[0], matTem);
	deMat(matTem);

	return RES_OK;
}



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
	BYTE *ptrOfMatRet, *ptrOfMatX, *ptrOfMatY;
	for (row = 0; row < matX->dim_row; ++row){
		for (col = 0; col < matY->dim_row; ++col)
		{
			int bytesOfRX = bytesOfRow(matX->dim_col);
			int bytesOfRR = bytesOfRow(matY->dim_row);

			int cntsVect = col / LENGTH;
			int offset = col % LENGTH;
			ptrOfMatRet = matRet->vect + row * bytesOfRR + cntsVect;

			int i;
			BYTE vectTem;
			for (i = 0; i < bytesOfRX; ++i){
				ptrOfMatX = matX->vect + row * bytesOfRX + i;
				ptrOfMatY = matY->vect + col * bytesOfRX + i;
				vectTem = (*ptrOfMatX) & (*ptrOfMatY);
				vectTem = sumFromVect(vectTem);
				vectTem = vectTem >> offset;
				(*ptrOfMatRet) ^= vectTem;
			}
		}
	}

	return matRet;
}


/* Add operation, as same as XOR */
Mat *addWithMask(
	const Mat *matX,
	const Mat *matY
	)
{
	if (matX == NULL || matY == NULL) return NULL;

	if (matX->dim_row != matY->dim_row ||
		matX->dim_col != matY->dim_col) return NULL;
	if (matX->vect == NULL || matY->vect == NULL) return NULL;

	/* Memory allocate */
	Mat *retMat = newMat(matX->dim_row, matX->dim_col, NULL, 0x00);
	if (retMat == NULL) return NULL;

	/* Calculate */
	int i;
	int bytesOfMat = matX->dim_row * bytesOfRow(matX->dim_col);
	for (i = 0; i < bytesOfMat; ++i){
		retMat->vect[i] = matX->vect[i] ^ matY->vect[i];
	}

	return retMat;
}



/* Real bitAnd operation  */
/* After operation, input matrices will all be deallocated */
Mat **bitAndWithMask(
	const Mat **matX,
	const Mat **matY
	)
{
	if (matX == NULL || matY == NULL) return NULL;

	Mat **matZ = (Mat **)malloc(MASKD * sizeof(Mat *));
	if (matZ == NULL) return NULL;
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
			if (i != 0 && j != 0) matTij[index] = bitAnd(matX[i], matY[j]);
			else if (i == 0 && j == 0){
				matTem = tenserProduct(matX[i], matY[j]);
				matTij[index] = multiply(matTem, matHat);
				deMat(matTem);
			}
			else if (i == 0){
				matTem = tenserProduct(matX[i], matY[j]);
				matTij[index] = multiply(matTem, matGrave);
				deMat(matTem);
			}
			else{/*  j == 0  */
				matTem = tenserProduct(matX[i], matY[j]);
				matTij[index] = multiply(matTem, matAcute);
				deMat(matTem);
			}
		}
	}

	BYTE *randMat = (BYTE *)malloc(matX[0]->dim_row * sizeof(BYTE));
	if (randMat == NULL) return NULL;

	/*  get matrix R  */
	for (i = 0; i < MASKD; ++i){
		
		/* get R(i,i) */
		index = i * MASKD + i;
		matR[index] = newMat(matX[0]->dim_row, matX[0]->dim_col, NULL, 0x00);
		memmove(matR[index]->vect, matTij[index]->vect, matX[0]->dim_row * sizeof(BYTE));
		deMat(matTij[index]);
		
		for (j = i + 1; j < MASKD; ++j){
					
			/* get R(i,j) through random generating */
			index = i * MASKD + j;
			matR[index] = newMat(matX[0]->dim_row, matX[0]->dim_col, NULL, 0x00);
			if (matR[index] == NULL) return NULL;

			for (k = 0; k < matTij[index]->dim_row; ++k){
				randMat[k] = (BYTE)rand();
			}
			
			memmove(matX[0]->vect, randMat, matX[0]->dim_row * sizeof(BYTE));

			/* get R(j,i)  */
			matTem = add(matR[index], matTij[index]);
			deMat(matTij[index]);

			index = j * MASKD + i;
			matR[index] = add(matTij[index], matTem);

			deMat(matTij[index]);
			deMat(matTem);

			/*  get R(0,j) */
			if (i == 0){
				index = i * MASKD + j;
				matTem = matR[index];
				matR[index] = multiply(matTem, matA);
				deMat(matTem);
			}
		}
	}

	for (i = 0; i < MASKD; ++i){
		matZ[i] = matR[i * MASKD];
		for (j = 1; j < MASKD; ++j){
			index = i * MASKD + j;
			matTem = matZ[i];
			matZ[i] = add(matTem, matR[index]);

			deMat(matTem);
			deMat(matR[index]);
		}
	}

	/* Deallocating*/
	free(randMat);
	free(matTij);
	free(matR);

	return matZ;
}

#else
/* Matrix multiplication */
Mat *multiply(
	const Mat *matX,
	const Mat *matY
	)
{
	if (matX == NULL || matY == NULL) return NULL;
	if (matX->dim_col != matY->dim_row) return NULL;
	if (matX->vect == NULL || matY->vect == NULL) return NULL;

	/* Memory allocation */
	Mat *matRet = newMat(matX->dim_row, matY->dim_col, NULL, 0x00);
	if (matRet == NULL) return NULL;

	/* Preparing a transposed matrix first */
	Mat *matYT = transpose(matY);

	/* Multipling */
	int row, col;
	BYTE *ptrOfMatRet, *ptrOfMatX, *ptrOfMatYT;
	for (row = 0; row < matX->dim_row; ++row){
		for (col = 0; col < matY->dim_col; ++col)
		{
			int bytesOfRX = bytesOfRow(matX->dim_col);
			int bytesOfRR = bytesOfRow(matY->dim_col);

			int cntsVect = col / LENGTH;
			int offset = col % LENGTH;
			ptrOfMatRet = matRet->vect + row * bytesOfRR + cntsVect;

			int i;
			BYTE vectTem;
			for (i = 0; i < bytesOfRX; ++i){
				ptrOfMatX = matX->vect + row * bytesOfRX + i;
				ptrOfMatYT = matYT->vect + col * bytesOfRX + i;
				vectTem = (*ptrOfMatX) & (*ptrOfMatYT);
				vectTem = sumFromVect(vectTem);
				vectTem = vectTem >> offset;
				(*ptrOfMatRet) ^= vectTem;
			}
		}
	}

	deMat(matYT);
	return matRet;
}

#endif


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
	Mat *retMat = newMat(matX->dim_row, matX->dim_col, NULL, 0x00);
	if (retMat == NULL) return NULL;

	/* Calculate */
	int i;
	int bytesOfMat = matX->dim_row * bytesOfRow(matX->dim_col);
	for (i = 0; i < bytesOfMat; ++i){
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
	for (i = 0; i < bytesOfMat; ++i){
		retMat->vect[i] = matX->vect[i] & matY->vect[i];
	}

	return retMat;
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

	BYTE *ptrOfSubMat;
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
	BYTE *ptrOfBigMat;
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
