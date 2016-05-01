//
//  LSbox.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#include "LSbox.h"

/*=========================================================*/
/*   MARK: Declaration     */
/*=========================================================*/
#if MASK

#if DIM_A
extern
BYTE matA[DIM_A], matInvA[DIM_A], matTransA[DIM_A];

#if !DIVIDE && DIM_L == 16
extern
BYTE matAs[L_SIZE], matInvAs[L_SIZE], matTransAs[L_SIZE];
#endif

static
BYTE  matT[L_SIZE] = { 0 };
#endif /* DIM_A */

#endif /* MASK */

static const
BYTE matL[L_SIZE] = MAT_LV;
static const
BYTE keyR[KEY_SIZE] = KEY_RV;
static const
BYTE rdConst[CONST_SIZE] = CONSTR;



/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/

/* 4-bit SBOX */

static
Res sbox4b(
BYTE *s4bRes,
const BYTE *mats4b,
const int *dims4b// [4, DIM_L]
)
#if MASK
{
	Res res = RES_OK;

	BYTE rvect[4][MASKD][2] = { 0 };
	BYTE imd[4][MASKD][2] = { 0 };
	BYTE product[MASKD][2] = { 0 };
	BYTE tem[MASKD][4][2] = { 0 };

	int btsVect, btsMat, i;
	const int dimsB[4] = { 1, DIM_L, 1, DIM_L };

	if (s4bRes == NULL) return RES_INVALID_POINTER;
	if (dims4b[0] != 4) return RES_INVALID_DIMENSION;

	btsMat = dims4b[0] * bytesOfRow(dims4b[1]);
	memset(s4bRes, 0, btsMat * MASKD * sizeof(BYTE));

	for (i = 0; i < MASKD; ++i){
		int j;
		for (j = 0; j < 4; ++j){
			rvect[j][i][0] = mats4b[i * btsMat + j * 2];
			rvect[j][i][1] = mats4b[i * btsMat + j * 2 + 1];
		}
	}

	// a = 1 x 2 + 0
	// b = 2 x a + 3
	// c = 1 x b + 2
	// d = b x c + 1	
	btsVect = MASKD * 2;

	res = bitAndWithMask((BYTE *)product, (const BYTE *)rvect + 1 * btsVect, (const BYTE *)rvect + 2 * btsVect, dimsB);
	CHECK(res);
	res = addWithMask((BYTE *)imd, (const BYTE *)product, (const BYTE *)rvect, dimsB);
	CHECK(res);

	res = bitAndWithMask((BYTE *)product, (const BYTE *)rvect + 2 * btsVect, (const BYTE *)imd, dimsB);
	CHECK(res);
	res = addWithMask((BYTE *)imd + 1 * btsVect, (const BYTE *)product, (const BYTE *)rvect + 3 * btsVect, dimsB);
	CHECK(res);

	res = bitAndWithMask((BYTE *)product, (const BYTE *)rvect + 1 * btsVect, (const BYTE *)imd + 1 * btsVect, dimsB);
	CHECK(res);
	res = addWithMask((BYTE *)imd + 2 * btsVect, (const BYTE *)product, (const BYTE *)rvect + 2 * btsVect, dimsB);
	CHECK(res);

	res = bitAndWithMask((BYTE *)product, (const BYTE *)imd + 1 * btsVect, (const BYTE *)imd + 2 * btsVect, dimsB);
	CHECK(res);
	res = addWithMask((BYTE *)imd + 3 * btsVect, (const BYTE *)product, (const BYTE *)rvect + 1 * btsVect, dimsB);
	CHECK(res);

	/* Generate the correct order
	 * d a b c
	 * 3 0 1 2
	 */
	for (i = 0; i != MASKD; ++i){
		s4bRes[i * btsMat]     = imd[3][i][0];
		s4bRes[i * btsMat + 1] = imd[3][i][1];

		s4bRes[i * btsMat + 2] = imd[0][i][0];
		s4bRes[i * btsMat + 3] = imd[0][i][1];

		s4bRes[i * btsMat + 4] = imd[1][i][0];
		s4bRes[i * btsMat + 5] = imd[1][i][1];

		s4bRes[i * btsMat + 6] = imd[2][i][0];
		s4bRes[i * btsMat + 7] = imd[2][i][1];

	}
	return res;
}

#else /* Unmask */
{
	Res res = RES_OK;
	BYTE rvect[4][DIM_L / 8] = { 0 };
	BYTE imd[4][DIM_L / 8] = { 0 };
	BYTE product[DIM_L / 8] = { 0 };
	
	int j;
	int dimsB[4] = { 1, DIM_L, 1, DIM_L};
	int btsVect;

	if (dims4b[0] != 4) return RES_INVALID_DIMENSION;

	for (j = 0; j != 4; ++j){
		rvect[j][0] = mats4b[j * 2];
		rvect[j][1] = mats4b[j * 2 + 1];
	}

	/* a = 1 x 2 + 0
	 * b = 2 x a + 3
	 * c = 1 x b + 2
	 * d = b x c + 1	
	 */
	btsVect = 2;
	res = bitAnd((BYTE *)product, (const BYTE *)rvect + btsVect, (const BYTE *)rvect + 2 * btsVect, dimsB);
	CHECK(res);
	res = add((BYTE *)imd, (const BYTE *)product, (const BYTE *)rvect, dimsB);
	CHECK(res);

	res = bitAnd((BYTE *)product, (const BYTE *)rvect + 2 * btsVect, (const BYTE *)imd, dimsB);
	CHECK(res);
	res = add((BYTE *)imd + btsVect, (const BYTE *)product, (const BYTE *)rvect + 3 * btsVect, dimsB);
	CHECK(res);

	res = bitAnd((BYTE *)product, (const BYTE *)rvect + btsVect, (const BYTE *)imd + btsVect, dimsB);
	CHECK(res);
	res = add((BYTE *)imd + 2 * btsVect, (const BYTE *)product, (const BYTE *)rvect + 2 * btsVect, dimsB);
	CHECK(res);

	res = bitAnd((BYTE *)product, (const BYTE *)imd + btsVect, (const BYTE *)imd + 2 * btsVect, dimsB);
	CHECK(res);
	res = add((BYTE *)imd + 3 * btsVect, (const BYTE *)product, (const BYTE *)rvect + btsVect, dimsB);
	CHECK(res);


	/* Generate the correct order
	 * d a b c  // 3 0 1 2 
	 */
	s4bRes[0] = imd[3][0];
	s4bRes[1] = imd[3][1];

	s4bRes[2] = imd[0][0];
	s4bRes[3] = imd[0][1];

	s4bRes[4] = imd[1][0];
	s4bRes[5] = imd[1][1];

	s4bRes[6] = imd[2][0];
	s4bRes[7] = imd[2][1];

	return res;
}
#endif /* MASK */



/* LBOX */

static
Res  lboxes(
BYTE *matsLin
)
#if MASK
{
	Res res = RES_OK;
	int i;
	/* z[0] = T  *  x[0] */
	BYTE matTem[DIM_L] = { 0 };
	const int dimsM[4] = { DIM_S, DIM_L, DIM_L, DIM_L };
	int btsOfMat;

	btsOfMat = DIM_S * bytesOfRow(DIM_L);
	memcpy((BYTE *)matTem, (const BYTE*)matsLin, btsOfMat * sizeof(BYTE));
#if DIM_A
	res = multiply(matsLin, (const BYTE *)matTem, (const BYTE *)matT, dimsM);
#else /* DIM_A == 0 */
	res = multiply(matsLin, (const BYTE *)matTem, (const BYTE *)matL, dimsM);
#endif
	CHECK(res);
	for (i = 1; i != MASKD; ++i){
		memcpy((BYTE *)matTem, (const BYTE*)matsLin + i * btsOfMat, btsOfMat * sizeof(BYTE));
		res = multiply((BYTE *)matsLin + i * btsOfMat, (const BYTE *)matTem, (const BYTE *)matL, dimsM);
		CHECK(res);
	}
	return res;
}
#else /* Unmask  */
{
	Res res = RES_OK;
	const int dims[4] = { DIM_S, DIM_L, DIM_L, DIM_L };
	BYTE tem[DIM_L] = { 0 };
	memcpy((BYTE *)tem, matsLin, DIM_L);
	res = multiply(matsLin, tem, (const BYTE *)matL, dims);
	return res;
}
#endif /* MASK */

/* SBOX*/

/* DIM_S-bit S-box (masked)*/
static
Res  sboxes(
BYTE *matsSin
)
#if MASK
{
	int i;
	Res res = RES_OK;
	BYTE s4bRes[MASKD][KEY_SIZE] = { 0 };
	BYTE left[MASKD][KEY_SIZE] = { 0 };
	BYTE right[MASKD][KEY_SIZE] = { 0 };
	BYTE sum[MASKD][KEY_SIZE] = { 0 };

	BYTE *ptrOfL, *ptrOfR;
	int  btsHalf;
	int  theLast;
	const int dimsF[4] = { DIM_S / 2, DIM_L, DIM_S / 2, DIM_L };

	btsHalf = DIM_S / 2 * bytesOfRow(DIM_L);
	theLast = MASKD - 1;

	for (i = 0; i != MASKD; ++i){
		memcpy((BYTE *)left  + i * btsHalf, (const BYTE *)matsSin + (i * 2)     * btsHalf, btsHalf);
		memcpy((BYTE *)right + i * btsHalf, (const BYTE *)matsSin + (i * 2 + 1) * btsHalf, btsHalf);
	}

	/* Combine a bigger sbox, from 4-bit to 8-bit */
	ptrOfL = (BYTE *)left, ptrOfR = (BYTE *)right;
	/* Feistel Struct Begins */
	for (i = 0; i != FEISTEL; ++i)
	{
		BYTE tem[KEY_SIZE] = { 0 };
		memcpy((BYTE *)tem, (const BYTE*)ptrOfL + theLast * btsHalf, btsHalf * sizeof(BYTE));
		/* add key_r to the last masked component */
		res = add(ptrOfL + theLast * btsHalf, (const BYTE *)tem, (const BYTE *)keyR, dimsF);
		CHECK(res);
		/* pass a 4-bit sbox */
		res = sbox4b((BYTE *)s4bRes, (const BYTE *)ptrOfL, dimsF);
		CHECK(res);
		
		/* recover the left matrices */
		memcpy((BYTE *)ptrOfL + theLast * btsHalf, (const BYTE*)tem, btsHalf * sizeof(BYTE));

		/* do 'XOR' with the right matrix */
		res = addWithMask((BYTE *)sum, (const BYTE *)ptrOfR, (const BYTE *)s4bRes, dimsF);
		CHECK(res);
		/* exchange  */
		ptrOfR = ptrOfL;
		ptrOfL = (BYTE *)sum;
	}

	/* Catenate those vectors to a matrix to return */
	for (i = 0; i != MASKD; ++i){
		memcpy((BYTE *)matsSin + (i * 2)     * btsHalf, (const BYTE *)ptrOfL + i * btsHalf, btsHalf);
		memcpy((BYTE *)matsSin + (i * 2 + 1) * btsHalf, (const BYTE *)ptrOfR + i * btsHalf, btsHalf);
	}
	return res;
}

#else /* Unmask */
{
	/* Split matrix to rowsUpper and rowsLower */
	Res res = RES_OK;
	BYTE s4bRes[KEY_SIZE] = { 0 };
	BYTE leftT[KEY_SIZE] = { 0 };
	BYTE rightT[KEY_SIZE] = { 0 };
	BYTE sumT[KEY_SIZE] = { 0 };
	BYTE *sum, *left, *right;

	const int dimsF[4] = { DIM_S / 2, DIM_L, DIM_S / 2, DIM_L };
	int btsHalf;
	int i;

	sum = (BYTE *)sumT; left = (BYTE *)leftT; right = (BYTE *)rightT;

	btsHalf = DIM_S / 2 * bytesOfRow(DIM_L);
	memcpy(left , matsSin , btsHalf);
	memcpy(right, matsSin + btsHalf, btsHalf);
	
	
	/* Combine a bigger sbox, from 4-bit to 8-bit */
	/* Feistel Struct */
	for (i = 0; i != FEISTEL; ++i)
	{
		/* do 'XOR' with key_r */
		res = add(sum, (const BYTE *)left, (const BYTE *)keyR, dimsF);
		CHECK(res);
		/* pass a 4-bit sbox */
		res = sbox4b((BYTE *)s4bRes, (const BYTE *)sum, dimsF);
		CHECK(res);
		/* do 'XOR' with the right matrix */
		res = add(sum, (const BYTE *)right, (const BYTE *)s4bRes, dimsF);
		CHECK(res);
		/* exchange each side */
		right = left;
		left = sum;
	}

	/* Catenate those vectors to a matrix to return */
	memcpy(matsSin, left, btsHalf);
	memcpy(matsSin + btsHalf, right , btsHalf);

	return res;
}
#endif /* MASK */







#if DIM_A /* !DIVIDE && DIM_A */
Res  getMatT()

{
	/* Get Matrix T  */
	BYTE matRight[L_SIZE] = { 0 };
	Res res = RES_OK;
	const int dims[4] = { DIM_L, DIM_L, DIM_L, DIM_L };
#if DIM_L == 16
	res = multiply(matRight, (const BYTE *)matTransAs, (const BYTE *)matL    , dims);
	res = multiply(matT	   , (const BYTE *)matInvAs  , (const BYTE *)matRight, dims);
#else
	res = multiply(matRight, (const BYTE *)matTransA, (const BYTE *)matL    , dims);
	res = multiply(matT    , (const BYTE *)matInvA  , (const BYTE *)matRight, dims);
#endif
	return res;
}
#endif






/*=========================================================*/
/*   MARK: Public   Functions      */
/*=========================================================*/

Res encrypto(
	BYTE *cipher,
	const BYTE *plain,
	const BYTE *key
	)
#if MASK
{
	Res res = RES_OK;
	BYTE matMasked[MASKD][DIM_L] = { 0 };
	BYTE tem[DIM_L] = { 0 };
	int theLast;
	const int dims[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	int indexOfRound;
	int btsMat;
	if (plain == NULL || key == NULL) return RES_INVALID_POINTER;

	/* Encode PlainText */
	res = encode((BYTE *)matMasked, plain);
	CHECK(res);
	
	/* Add Key */
	theLast = MASKD - 1;
	btsMat = DIM_S * bytesOfRow(DIM_L);

	memcpy((BYTE *)tem, (const BYTE *)matMasked + theLast * btsMat, btsMat * sizeof(BYTE));
	res = add((BYTE *)matMasked + theLast * btsMat, key, (const BYTE *)tem, dims);
	
	for (indexOfRound = 0; indexOfRound != ROUNDS; ++indexOfRound)
	{
		BYTE roundK[DIM_L] = { 0 };

		res = sboxes((BYTE *)matMasked);
		CHECK(res);
		res = lboxes((BYTE *)matMasked);
		CHECK(res);
		res = add((BYTE *)roundK, rdConst, key, dims);
		CHECK(res);
		/* Add Key And Round Constant */
		memcpy((BYTE *)tem, (const BYTE *)matMasked + theLast * btsMat, DIM_L);
		res = add((BYTE *)matMasked + theLast * btsMat, (const BYTE *)tem, (const BYTE *)roundK, dims);
		CHECK(res);
	}
	/* Decode Cipher */
	res = decode(cipher, (const BYTE *)matMasked);
	return res;

}
#else
{
	Res  res = RES_OK;
	BYTE sum[DIM_L] = { 0 };
	
	int btsMat, round_i;
	const int dims[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	
	if (plain == NULL || key == NULL || ROUNDS < 0) return RES_INVALID_POINTER;
	btsMat = DIM_S * bytesOfRow(DIM_L);
	res = add(cipher, plain, key, dims);
	CHECK(res);

	for (round_i = 0; round_i != ROUNDS; ++round_i)
	{
		res = sboxes(cipher);
		CHECK(res);
		res = lboxes(cipher);
		CHECK(res);

		res = add(sum, cipher, key, dims);
		CHECK(res);
		res = add(cipher, (const BYTE *)sum, rdConst, dims);
		CHECK(res);
	}

	return res;
}
#endif /* MASK */

Res encrypto_fixed(){
	Res res = RES_OK;
	BYTE cipher[DIM_L] = { 0 };

	const BYTE plainT[DIM_L] = { 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf, 0xa0 };
	const BYTE cipherK[DIM_L] = { 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf, 0xb0 };

	res = encrypto(cipher, plainT, cipherK);
	return res;
}