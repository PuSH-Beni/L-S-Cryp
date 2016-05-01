//
//  Fundamentals.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//


#include "Fundamentals.h"

/*=========================================================*/
/*   MARK: Declaration     */
/*=========================================================*/
#if DIM_A
/* DIM_A == 4 or DIM_A == 8 */
static
BYTE matHat[DIM_A * DIM_A] = { 0 }, matGrave[DIM_A * DIM_A] = { 0 }, matAcute[DIM_A * DIM_A] = { 0 };

BYTE matA[DIM_A] = { 0 }, matInvA[DIM_A] = { 0 }, matTransA[DIM_A] = { 0 };
BYTE matAs[L_SIZE] = { 0 }, matInvAs[L_SIZE] = { 0 }, matTransAs[L_SIZE] = { 0 };
#endif


/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/


/* Calculate the # of bytes in one row */
int bytesOfRow(
int col
)
{
	int bytes;
	// bytes of each row :: if dim_col < LENGTH, allocate a byte as well
	if (col < LENGTH) return 1;

	#if LENG16
	bytes = col >> 4;
	#elif LENG8
	bytes = col / 8; // col / LENGTH
	#endif
	bytes = (col % LENGTH) ? bytes + 1 : bytes;
	return bytes;
}




/* Shift bit form j -> i */
static
BYTE shiftBit(
BYTE orig,
int i,
int j
)
{
	BYTE tem;
	tem = (UNIT_BYTE >> j) & orig;

	if (j > i) tem <<= (j - i);
	else if (i > j) tem >>= (i - j);
	else return tem;

	return tem;
}



/* Get sum of a byte */

static
BYTE sumOfByte(
BYTE byte
)
{
	BYTE ret;
	ret = 0;
	while (byte > 0){
		ret ^= UNIT_BYTE;
		byte &= (byte - 1);
	}

	return ret;
}



/* Transposition */
//static
Res transpose(
BYTE *transRes,
const BYTE *matOrig,
const int *dims
)
{
	int colOrig, rowOrig;
	int cntBytesOrig, cntBytesRet;
	int offsetOrig, offsetRet;
	BYTE vectTransed;
	BYTE *ptrOfVectRet;
	const BYTE *ptrOfVectOrig;

	int bytesOfRowOrig;
	int bytesOfRowRet;

	bytesOfRowOrig = bytesOfRow(dims[1]);
	bytesOfRowRet = bytesOfRow(dims[0]);
	memset(transRes, 0, dims[0] * bytesOfRowOrig);
	if (matOrig == NULL || transRes == NULL) return RES_INVALID_POINTER;
	/* Transposing */
	for (colOrig = 0; colOrig < dims[1]; ++colOrig){
		cntBytesOrig = colOrig / LENGTH;
		offsetOrig = colOrig % LENGTH;
		for (rowOrig = 0; rowOrig < dims[0]; ++rowOrig){
			/* the bit (rowOrig, colOrig) */
			cntBytesRet = rowOrig / LENGTH;
			offsetRet = rowOrig % LENGTH;

			ptrOfVectOrig = matOrig + bytesOfRowOrig * rowOrig + cntBytesOrig;
			ptrOfVectRet = transRes + bytesOfRowRet * colOrig + cntBytesRet;

			vectTransed = shiftBit(*ptrOfVectOrig, offsetRet, offsetOrig);
			(*ptrOfVectRet) ^= vectTransed;
		}
	}
	return RES_OK;

}



/* Generate DIM_A x DIM_A random matrices A, inv(A), trans(A)  */
#if DIM_A
static
Res genRandMat(
)
{
	Res res;
	int i;
	const int dimsE[4] = { DIM_A, DIM_A, DIM_A, DIM_A };
	BYTE matP[DIM_A] = { 0 };
	BYTE matE[DIM_A] = UNIT_MAT;
	BYTE tem[DIM_A] = { 0 };

	/*  get matrix A and inv(A) temporarily */
	memcpy((BYTE *)matTransA, (const BYTE *)matE, DIM_A * sizeof(BYTE));
	memcpy((BYTE *)matInvA, (const BYTE *)matE, DIM_A * sizeof(BYTE));

	/* get matrix P */

	for (i = 0; i < DIM_A; ++i){

		BYTE rowToAdd;
		BYTE zeroToSet;

		rowToAdd = (BYTE)rand();
		zeroToSet = UNIT_BYTE >> i;
		rowToAdd &= (~zeroToSet);
		#if DIM_A == 4
		rowToAdd &= 0xf0;
		#endif
		matE[i] ^= rowToAdd;

		/* Tanspose(P) */
		
		res = transpose(matP, (const BYTE *)matE, dimsE);
		CHECK(res);
		/*  A = P x A   ==>  A^T = A^T x P^T  */
		
		memcpy((BYTE *)tem, (const BYTE *)matTransA, DIM_A * sizeof(BYTE));
		res = multiply(matTransA, (const BYTE *)tem, (const BYTE *)matE, dimsE);
		CHECK(res);
		/*  A^{-1} = A^{-1} x P     */
		memcpy((BYTE *)tem, (const BYTE *)matInvA, DIM_A * sizeof(BYTE));
		res = multiply(matInvA, (const BYTE *)tem, (const BYTE *)matP, dimsE);
		CHECK(res);
		/* Refresh matE to Unit Matrix */
		matE[i] ^= rowToAdd;
	}

	res = transpose(matA, (const BYTE *)matTransA, dimsE);
	return res;
}




/* Generate hatA, graveA and acuteA */
/* All of them are DIM_A x DIM_A^2 */
static
Res hatA(
)
{
	int i, j;
	Res res = RES_OK;
	BYTE *ptrOfMatATi;
	BYTE *ptrOfMatATj, *ptrOfMatR;
	const int dimsH[4] = { DIM_A, DIM_A, DIM_A * DIM_A, DIM_A };

	if (matTransA == NULL || matInvA == NULL) return RES_NON_MATA;
	//row: DIM_A^2, col: DIM_A
	BYTE matRight[DIM_A * DIM_A * DIM_A / LENGTH] = { 0 };
	/* get transposition of the right matrix: (E' x (A tp B))^T */
	for (i = 0; i < DIM_A; ++i){
		ptrOfMatATi = matTransA + i;
		for (j = 0; j < DIM_A; ++j){
			ptrOfMatATj = matTransA + j;
			ptrOfMatR = matRight + i * DIM_A + j;
			(*ptrOfMatR) = (*ptrOfMatATi) & (*ptrOfMatATj);
		}
	}

	/* get \hat{A} now */	
	res = multiply(matHat, (const BYTE *)matInvA, (const BYTE *)matRight, dimsH);
	return res;
}



static
Res acuteA(
)
{
	int i, j;
	Res res = RES_OK;
	BYTE matRight[DIM_A * DIM_A * DIM_A / LENGTH] = { 0 };
	BYTE *ptrOfMatATj, *ptrOfMatR, *ptrOfMatI;
	BYTE matE[DIM_A * DIM_A / LENGTH] = UNIT_MAT;
	const int dimsH[4] = { DIM_A, DIM_A, DIM_A * DIM_A, DIM_A };
	if (matTransA == NULL || matInvA == NULL) return RES_NON_MATA;

	/* get transposition of the right matrix: (E' x (E tp A))^T */
	for (i = 0; i < DIM_A; ++i){
		ptrOfMatI = matE + i;
		for (j = 0; j < DIM_A; ++j){
			ptrOfMatATj = matTransA + j;
			ptrOfMatR = matRight + i * DIM_A + j;
			(*ptrOfMatR) = (*ptrOfMatI) & (*ptrOfMatATj);
		}
	}

	/* get \grave{A} now */
	res = multiply(matAcute, (const BYTE *)matInvA, (const BYTE *)matRight, dimsH);
	return res;
}



static
Res graveA(
)
{
	int i, j;
	Res res = RES_OK;
	BYTE *ptrOfMatATi, *ptrOfMatR, *ptrOfMatI;
	BYTE matE[DIM_A * DIM_A / LENGTH] = UNIT_MAT;

	BYTE matRight[DIM_A * DIM_A * DIM_A / LENGTH] = { 0 };
	const int dimsH[4] = { DIM_A, DIM_A, DIM_A * DIM_A, DIM_A };
	if (matTransA == NULL || matInvA == NULL) return RES_NON_MATA;
	/* get transposition of the right matrix: (E' x (A tp E))^T */
	for (i = 0; i < DIM_A; ++i){
		ptrOfMatATi = matTransA + i;
		for (j = 0; j < DIM_A; ++j){
			ptrOfMatI = matE + j;
			ptrOfMatR = matRight + i * DIM_A + j;
			(*ptrOfMatR) = (*ptrOfMatI) & (*ptrOfMatATi);
		}
	}

	/* get \Acute{A} now */	
	res = multiply(matGrave, (const BYTE *)matInvA, (const BYTE *)matRight, dimsH);
	return res;
}

#endif /* DIM_A != 0 */


/* Tensor Product(kron) for two vectors 
 * Return a n x n^2 matrix
 * row: 1, col: X_COL^2;  BYTE[DIM_S*X_ROW]
 */
static
Res tensorProduct(
BYTE *tensProdRes,
const BYTE *matX,
const BYTE *matY,
const int *dims // {1, DIM_A}
)
{
	int i, j;

	if (matX == NULL || matY == NULL || tensProdRes == NULL) return RES_INVALID_POINTER;
	if ((dims[0] != dims[2]) ||
	(dims[1] != dims[3]) ||
	(dims[1] > 8))return RES_INVALID_DIMENSION;

	memset(tensProdRes, 0, dims[0] * bytesOfRow(dims[1]));

	#if DIM_A == 4
	int bts = bytesOfRow(X_COL * X_COL);
	for (i = 0; i < X_ROW; ++i) {
		BYTE ident = UNIT_BYTE;

		for (j = 0; j < X_COL ; j += 2){
			BYTE tem = 0x00;
			if ((BYTE)matX[i] & ident){
				tem ^= matY[i];
			}
			ident = ident >> 1;
			if ((BYTE)matX[i] & ident){
				tem ^= (BYTE)matY[i] >> 4;
			}
			matR[i * bts + j / 2] = tem;
			ident >>= 1;

		}
	}
	#else /* DIM_A == 8 or 0 */

	for (i = 0; i != dims[0]; ++i){
		BYTE unit = UNIT_BYTE;
		for (j = 0; j != dims[3]; ++j){
			tensProdRes[i * dims[1] + j] = (matX[i] & unit) ? matY[i] : 0x00;
			unit >>= 1;
		}
	}

	#endif /* DIM_A == 4 */
	return RES_OK;
}



/*=========================================================*/
/*   MARK: Public   Functions      */
/*=========================================================*/
#if MASK
#if DIM_A
/* Set-up */
Res setupEnc(
)
{
	int btsOfRow, btsOfMat;
	int indexOfSlice, indexOfByte, indexOfDest;
	/* Get Matrix A, inv(A), trans(A) */
	Res res;


	res = genRandMat();
	CHECK(res);
	/* Get Matrix \hat{A}, \grave{A}, \acute{A} */
	res = hatA();
	CHECK(res);
	res = graveA();
	CHECK(res);
	res = acuteA();
	CHECK(res);

	#if !DIVIDE
	/*| inv(A), 0 |
	| 0, inv(A) | */

	btsOfRow = bytesOfRow(DIM_L);
	btsOfMat = DIM_A * bytesOfRow(DIM_A);
	indexOfDest = 0;
	for (indexOfSlice = 0; indexOfSlice != DIVIDE_PARTS; ++indexOfSlice){
		#if DIM_A == 4
		BYTE oddFlag = (BYTE)indexOfSlice & 0x01;
		#endif
		for (indexOfByte = 0; indexOfByte != btsOfMat; ++indexOfByte){
			#if DIM_A == 4
			matAs[indexOfDest]      = oddFlag ? matA[indexOfByte] >> 4
			: matA[indexOfByte];
			matInvAs[indexOfDest]   = oddFlag ? matInvA[indexOfByte] >> 4
			: matInvA[indexOfByte];
			matTransAs[indexOfDest] = oddFlag ? matTransA[indexOfByte] >> 4
			: matTransA[indexOfByte];
			#else
			matAs[indexOfDest] = matA[indexOfByte];
			matInvAs[indexOfDest] = matInvA[indexOfByte];
			matTransAs[indexOfDest] = matTransA[indexOfByte];
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

	return res;
}

#endif /* DIM_A */

/* Encode(Mask) the plain text */
Res encode(
BYTE *matMasked,
const BYTE *matPlain
)
{
	int i, j;
	Res res;
	int btsOfMat;
	BYTE randMat[DIM_L] = { 0 };
	BYTE matSumOfRest[DIM_L] = { 0 };
	BYTE tem[DIM_L] = { 0 };
	const int dimsM[4] = { DIM_S, DIM_L, DIM_L, DIM_L };
	const int dimsAdd[4] = { DIM_S, DIM_L, DIM_S, DIM_L };

	if (matPlain == NULL || matMasked == NULL) return RES_INVALID_POINTER;
	btsOfMat = DIM_S *  bytesOfRow(DIM_L);
	memset(matMasked, 0, btsOfMat * MASKD);

	for (i = 1; i < MASKD; ++i){
		for (j = 0; j < btsOfMat; ++j){
			randMat[j] = (BYTE)rand();
		}
		if (matMasked + i * btsOfMat == NULL) return RES_INVALID_POINTER;
		memcpy(matMasked + i * btsOfMat, (const BYTE *)randMat, btsOfMat * sizeof(BYTE));

		if (i == 1){
			memcpy((BYTE *)matSumOfRest, (const BYTE *)matMasked + i * btsOfMat, btsOfMat * sizeof(BYTE));
		}
		else {
			memcpy(tem, (const BYTE *)matSumOfRest, btsOfMat * sizeof(BYTE));
			res = add((BYTE *)matSumOfRest, (const BYTE *)tem, (const BYTE *)matMasked + i * btsOfMat, dimsAdd);
			CHECK(res);
		}

	}

	#if DIM_A

	#if DIM_L == 16 && !DIVIDE
	res = add(tem, (const BYTE *)matSumOfRest, matPlain, dimsAdd);
	CHECK(res);
	/*  = (x^T) x inv(A) ^ T */
	res = multiply(matMasked, (const BYTE *)tem, (const BYTE *)matInvAs, dimsM);
	#else
	res = add(tem, matSumOfRest, matPlain, dimsAdd);
	/*  = (x^T) x inv(A) ^ T */
	const int dims[4] = {DIM_S, DIM_L, DIM_A, DIM_A};
	res = multiply(matMasked[0], (const BYTE *)tem, (const BYTE *)matInvA, dims);

	#endif

	#else /* DIM_A */
	res = add(matMasked, (const BYTE *)matSumOfRest, matPlain, dimsAdd);
	#endif
	return res;
}


/* Decode(Unmask) the Secret text */
Res decode(
BYTE *matUnmask,
const BYTE *matsSecret
)
{
	int i;
	Res res = RES_OK;
	const int dimsM[4] = { DIM_S, DIM_L, DIM_L, DIM_L };
	const int dimsAdd[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	if (matsSecret == NULL || matUnmask == NULL) return RES_INVALID_POINTER;
	memset(matUnmask, 0, DIM_L * sizeof(BYTE));
	#if DIM_A

	#if !DIVIDE && DIM_L == 16
	res = multiply(matUnmask, matsSecret, (const BYTE *)matAs, dimsM);
	CHECK(res);
	#else // DIM_L != 16
	const int dims[4] = {DIM_S, DIM_L, DIM_A, DIM_A};
	res = multiply(matUnmask, matsSecret, (const BYTE *)matA, dims);
	#endif

	#else /* DIM_A == 0 */
	memcpy(matUnmask, matsSecret, DIM_L * sizeof(BYTE));
	#endif
	for (i = 1; i < MASKD; ++i){
		BYTE tem[DIM_L] = { 0 };
		memcpy((BYTE *)tem, (const BYTE *)matUnmask, DIM_L * sizeof(BYTE));
		res = add(matUnmask, (const BYTE *)tem, matsSecret + i * DIM_L, dimsAdd);
		CHECK(res);
	}
	return res;
}


/* Add operation, as same as XOR */
Res addWithMask(
BYTE *addRes,
const BYTE *matEX,
const BYTE *matEY,
const int *dims
)
{
	int i;
	Res res = RES_OK;
	int bytsMat;

	if (matEX == NULL || matEY == NULL || addRes == NULL) return RES_INVALID_POINTER;
	if (dims[0] != dims[2] || dims[1] != dims[3]) return RES_INVALID_DIMENSION;
	/* Calculate */
	bytsMat = dims[0] * bytesOfRow(dims[1]);
	memset(addRes, 0, bytsMat * MASKD);
	for (i = 0; i < MASKD; ++i){
		res = add(addRes + i * bytsMat, matEX + i * bytsMat, matEY + i * bytsMat, dims);
		CHECK(res);
	}
	return res;
}



/* Masked bitAnd operation  */
Res bitAndWithMask(
BYTE *bitAndRes, // bitAndRes[MASKD][2]
const BYTE *matEX,
const BYTE *matEY,
const int *dims //{1, DIM_L, 1, DIM_L}
)
#if DIM_L == 16 && DIM_A
{
	int i, j;
	int index;

	int slicesIndex;
	Res res = RES_OK;
	BYTE matTij[MASKD_SQURE] = { 0 };
	BYTE matR[MASKD_SQURE] = { 0 };

	BYTE matEZPart[MASKD] = { 0 };
	BYTE matEYPart[MASKD] = { 0 };
	BYTE matEXPart[MASKD] = { 0 };

	BYTE matTem[DIM_A] = { 0 };
	int btsMat;

	if (matEX == NULL || matEY == NULL || bitAndRes == NULL) return RES_INVALID_POINTER;
	if (dims[0] != 1 || dims[0] != dims[2] || dims[1] != dims[3]) return RES_INVALID_DIMENSION;
	btsMat = dims[0] * bytesOfRow(dims[1]);
	memset(bitAndRes, 0, btsMat * MASKD);

	for (slicesIndex = 0; slicesIndex != DIVIDE_PARTS; ++slicesIndex){
		#if DIM_A == 4
		BYTE oddFlag = (BYTE)slicesIndex & 0x01;
		#endif
		for (i = 0; i != MASKD; ++i){
			#if DIM_A == 8
			matEXPart[i] = matEX[i * btsMat + slicesIndex];
			matEYPart[i] = matEY[i * btsMat + slicesIndex];
			#elif DIM_A == 4
			int offset = (slicesIndex < 2) ? 0 : 1;
			matEXPart[i] = matEX[i][offset];
			matEYPart[i] = matEY[i][offset];

			if (oddFlag){//odd
				matEXPart[i] &= 0x0f;
				matEXPart[i] <<= 4;

				matEYPart[i] &= 0x0f;
				matEYPart[i] <<= 4;
			}
			else{
				matEXPart[i] &= 0xf0;
				matEYPart[i] &= 0xf0;
			}
			#endif
		}
		/*  get matrix T  */
		for (i = 0; i != MASKD; ++i){
			for (j = 0; j != MASKD; ++j){
				const int dimsM[4] = { 1, DIM_A*DIM_A, DIM_A, DIM_A*DIM_A };
				const int dimsTensor[4] = { 1, DIM_A, 1, DIM_A };
				index = i * MASKD + j;
				#if DIM_A
				if (i == 0 && j == 0){
					res = tensorProduct(matTem, (const BYTE *)matEXPart + i, (const BYTE *)matEYPart + j, dimsTensor);
					CHECK(res);
					res = multiply(matTij + index, (const BYTE *)matTem, (const BYTE *)matHat, dimsM);
					CHECK(res);
				}
				else if (i == 0){
					res = tensorProduct(matTem, (const BYTE *)matEXPart + i, (const BYTE *)matEYPart + j, dimsTensor);
					CHECK(res);
					res = multiply(matTij + index, (const BYTE *)matTem, (const BYTE *)matGrave, dimsM);
					CHECK(res);
				}
				else if (j == 0){
					res = tensorProduct(matTem, (const BYTE *)matEXPart + i, (const BYTE *)matEYPart + j, dimsTensor);
					CHECK(res);
					res = multiply(matTij + index, (const BYTE *)matTem, (const BYTE *)matAcute, dimsM);
					CHECK(res);
				}
				else {/*  i != 0 && j != 0  */
					res = bitAnd(matTij + index, (const BYTE *)matEXPart + i, (const BYTE *)matEYPart + j, dimsTensor);
					CHECK(res);
				}

				#else  /* DIM_A == 0 */
				res = bitAnd(matTij + index, (const BYTE *)matEXPart + i, (const BYTE *)matEYPart + j, dimsTensor);
				#endif /* DIM_A */
			}
		}
		/*  get matrix R  */
		for (i = 0; i < MASKD; ++i){
			/* get R(i,i) */
			index = i * MASKD + i;
			matR[index] = matTij[index];
			for (j = i + 1; j < MASKD; ++j){
				BYTE tem;
				int indexJI;
				/* get R(i,j) through random generating */
				index = i * MASKD + j;
				indexJI = j * MASKD + i;
				#if DIM_A == 4
				randMat &= 0xf0;
				#endif
				matR[index] = (BYTE)rand();
				/* get R(j,i)  */	
				matR[indexJI] = matTij[indexJI] ^ matR[index] ^ matTij[index];
				#if DIM_A
				/*  get R(0,j) */
				if (i == 0){
					index = j;
					tem = matR[index];
					const int dimsT[4] = { 1, DIM_A, DIM_A, DIM_A };
					res = multiply(matR + index, (const BYTE *)&tem, (const BYTE *)matA, dimsT);
					CHECK(res);
				}
				#endif  /* DIM_A != 0 */
			}
		}

		for (i = 0; i != MASKD; ++i){
			matEZPart[i] = matR[i];
			for (j = 1; j != MASKD; ++j){
				index = j * MASKD + i;
				matEZPart[i] ^= matR[index];
			}
		}
		for (i = 0; i != MASKD; ++i){
			#if DIM_A == 8
			bitAndRes[i * btsMat + slicesIndex] = matEZPart[i];
			#elif DIM_A == 4
			int offset = (slicesIndex < 2) ? 0 : 1;
			if (oddFlag){//odd
				BYTE tem = matEZPart[i] >> 4;
				bitAndRes[i][offset] ^= (tem & 0x0f);
			}
			else bitAndRes[i][offset] ^= matEZPart[i];

			#endif
		}
	}
	return res;
}

#else
{
	Res res = RES_OK;
	BYTE matTij[MASKD_SQURE][DIM_L / 8] = { 0 };
	BYTE matR[MASKD_SQURE][DIM_L / 8] = { 0 };
	int i, j, k;
	int index, btsMat;

	if (matEX == NULL || matEY == NULL) return RES_INVALID_POINTER;
	btsMat = dims[0] * bytesOfRow(dims[1]);
	if (btsMat != DIM_L / 8) return RES_INVALID_DIMENSION;

	/*  get matrix T  */
	for (i = 0; i != MASKD; ++i){
		for (j = 0; j != MASKD; ++j){
			index = i * MASKD + j;
			#if DIM_A
			if (i == 0 && j == 0){
				res = tensorProduct(matTem, (const BYTE *)matEX + i * btsMat, (const BYTE *)matEY + j * btsMat, dims);
				CHECK(res);
				res = multiply(matTij + index * btsMat, (const BYTE *)matTem, (const BYTE *)matHat, dimsM);
				CHECK(res);
			}
			else if (i == 0){
				res = tensorProduct(matTem, (const BYTE *)matEX+i * btsMat, (const BYTE *)matEY+j * btsMat, dims);
				CHECK(res);
				res = multiply(matTij + index * btsMat, (const BYTE *)matTem, (const BYTE *)matGrave, dimsM);
				CHECK(res);
			}
			else if (j == 0){
				res = tensorProduct(matTem, (const BYTE *)matEX+i* btsMat, (const BYTE *)matEY+j * btsMat, dims);
				CHECK(res);
				res = multiply(matTij + index * btsMat, (const BYTE *)matTem, (const BYTE *)matAcute, dimsM);
				CHECK(res);
			}
			else {/*  i != 0 && j != 0  */
				res = bitAnd(matTij + index, (const BYTE *)matEX + i, (const BYTE *)matEY + j, dims);
				CHECK(res);
		}

			#else  /* DIM_A == 0 */
			//res = bitAnd((BYTE *)matTij + index * btsMat, (const BYTE *)matEX + i * btsMat, (const BYTE *)matEY + j * btsMat, dims);
			//CHECK(res);
			matTij[index][0] = matEX[i * btsMat]	 & matEY[j * btsMat];
			matTij[index][1] = matEX[i * btsMat + 1] & matEY[j * btsMat + 1];
#endif /* DIM_A */
		}
	}

	/*  get matrix R  */
	for (i = 0; i != MASKD; ++i){
		/* get R(i,i) */
		index = i * MASKD + i;
		matR[index][0] = matTij[index][0];
		matR[index][1] =  matTij[index][1];
		for (j = i + 1; j != MASKD; ++j){
			BYTE tem[DIM_L / 8] = { 0 };
			/* get R(i,j) through random generating */
			index = i * MASKD + j;
#if DIM_A == 4
				randMat[0] &= 0xf0;
				randMat[1] &= 0xf0;
#endif
			matR[index][0] = (BYTE)rand();
			matR[index][1] = (BYTE)rand();

			/* get R(j,i)  */
			tem[0] = matR[index][0] ^ matTij[index][0];
			tem[1] = matR[index][1] ^ matTij[index][1];
			
			index = j * MASKD + i;
			matR[index][0] = matTij[index][0] ^ tem[0];
			matR[index][1] = matTij[index][1] ^ tem[1];

			#if DIM_A
			/*  get R(0,j) */
			if (i == 0){
				index = j;
				matTem = matR[index];
				matR[index] = multiply(matTem, matA);
				deMat(matTem);
			}
			#endif  /* DIM_A */
		}
	}

	for (i = 0; i != MASKD; ++i){
		bitAndRes[i * btsMat]	  = matR[i][0];
		bitAndRes[i * btsMat + 1] = matR[i][1];
		for (j = 1; j != MASKD; ++j){
			index = j * MASKD + i;
			bitAndRes[i * btsMat]	  ^= matR[index][0];
			bitAndRes[i * btsMat + 1] ^= matR[index][1];
		}
	}
	return res;
}
#endif /* DIM_L == 16 && DIM_A */

#endif /* MASK */

/*  Matrix Multiply(Transposed), (A, B) == A x B^T   */
Res multiply(
BYTE *multiRes, //dims[0]*dims[2]
const BYTE *matX,
const BYTE *matY,
const int *dims
)
{
	int row, col;
	BYTE *ptrOfMatRet;
	const BYTE *ptrOfMatX, *ptrOfMatY;
	int bytesOfRX, bytesOfRR;
	
	if (matX == NULL || matY == NULL || multiRes == NULL) return RES_INVALID_POINTER;
	if (dims[1] != dims[3]) return RES_INVALID_DIMENSION;

	bytesOfRX = bytesOfRow(dims[1]);
	bytesOfRR = bytesOfRow(dims[2]);
	memset(multiRes, 0, dims[0] * bytesOfRR);
	/* Multiplying */
	for (row = 0; row < dims[0]; ++row){
		for (col = 0; col < dims[2]; ++col)
		{
			int i;
			BYTE vectTem;
			int cntsVect, offset;


			cntsVect = col / LENGTH;
			offset = col % LENGTH;
			ptrOfMatRet = multiRes + row * bytesOfRR + cntsVect;

			for (i = 0; i < bytesOfRX; ++i){
				ptrOfMatX = matX + row * bytesOfRX + i;
				ptrOfMatY = matY + col * bytesOfRX + i;

				//vectTem = (BYTE)__builtin_popcount((*ptrOfMatX) & (*ptrOfMatY));
				//vectTem <<= (LENGTH - 1 - offset);
				vectTem = sumOfByte((*ptrOfMatX) & (*ptrOfMatY));
				vectTem >>= offset;
				(*ptrOfMatRet) ^= vectTem;
			}
		}
	}

	return RES_OK;
}


/* Simple Add operation, as same as XOR */
Res add(
BYTE *addRes,
const BYTE *matX,
const  BYTE *matY,
const int *dims
)
{
	int i;
	int bytesOfMat;

	if (matX == NULL || matY == NULL || addRes == NULL) return RES_INVALID_POINTER;
	if (dims[0] != dims[2] ||
	dims[1] != dims[3]) return RES_INVALID_DIMENSION;

	bytesOfMat = dims[0] * bytesOfRow(dims[1]);
	memset(addRes, 0, bytesOfMat);
	/* Calculate */
	for (i = 0; i != bytesOfMat; ++i){
		addRes[i] = matX[i] ^ matY[i];
	}

	return RES_OK;
}



/* Simple bitAnd operation , as same as AND */
Res bitAnd(
BYTE *bitAndRes,
const BYTE *matX,
const BYTE *matY,
const int *dims
)
{
	int i;
	int bytesOfMat;


	if (matX == NULL || matY == NULL || bitAndRes == NULL) return RES_INVALID_POINTER;
	if (dims[0] != dims[2] ||
	dims[1] != dims[3]) return RES_INVALID_DIMENSION;

	bytesOfMat = dims[0] * bytesOfRow(dims[1]);
	memset(bitAndRes, 0, bytesOfMat);
	/* calculate */	
	for (i = 0; i != bytesOfMat; ++i){
		bitAndRes[i] = matX[i] & matY[i];
	}

	return RES_OK;
}




/* =========================================== */
/*                The   End                    */
/* =========================================== */
