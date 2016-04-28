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
#if DIVIDE
static 
BYTE **keyRoundSlices, **rdConstSlices[ROUNDS], **LSlices, *TSlices[DIVIDE_PARTS];
#endif

#if DIM_A
extern
BYTE matA[DIM_A], matInvA[DIM_A], matTransA[DIM_A];

#if !DIVIDE && DIM_L == 16
extern
BYTE matAs[DIM_L * 2], matInvAs[DIM_L * 2], matTransAs[DIM_L * 2];
#endif

static
BYTE  matT[DIM_L] = { 0 };
#endif /* DIM_A */

#endif /* MASK */
/********************************************
 *  static constants
 ********************************************/
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
#if MASK
static
Res sbox4b(
BYTE *s4bRes,
const BYTE *mats4b,
const int *dims4b// [4, DIM_L]
)
{
	int i;

	Res res = RES_OK;

	BYTE tem[MASKD][4][2] = { 0 };
	BYTE rvectWithMask[4][MASKD][2] = { 0 };
	BYTE imdWithMask[4][MASKD][2] = { 0 };
	BYTE product[MASKD][2] = { 0 };
	int btsVect;
	int btsMat;
	const int dimsB[4] = { 1, DIM_L, 1, DIM_L };

	
	btsMat = dims4b[0] * bytesOfRow(dims4b[1]);

	if (s4bRes == NULL) return RES_INVALID_POINTER;
	if (dims4b[0] != 4) return RES_INVALID_DIMENSION;
	for (i = 0; i < MASKD; ++i){
		int j;
		for (j = 0; j < 4; ++j){
			rvectWithMask[j][i][0] = mats4b[i * btsMat + j * 2];
			rvectWithMask[j][i][1] = mats4b[i * btsMat + j * 2 + 1];
		}
	}

	// a = 1 x 2 + 0
	// b = 2 x a + 3
	// c = 1 x b + 2
	// d = b x c + 1
	
	
	btsVect = MASKD * 2;
	res = bitAndWithMask((BYTE *)product, (const BYTE *)rvectWithMask + 1 * btsVect, (const BYTE *)rvectWithMask + 2 * btsVect, dimsB);
	res = addWithMask((BYTE *)imdWithMask, (const BYTE *)product, (const BYTE *)rvectWithMask, dimsB);

	res = bitAndWithMask((BYTE *)product, (const BYTE *)rvectWithMask + 2 * btsVect, (const BYTE *)imdWithMask, dimsB);
	res = addWithMask((BYTE *)imdWithMask+ 1 * btsVect, (const BYTE *)product, (const BYTE *)rvectWithMask + 3 * btsVect, dimsB);

	res = bitAndWithMask((BYTE *)product, (const BYTE *)rvectWithMask + 1 * btsVect, (const BYTE *)imdWithMask + 1 * btsVect, dimsB);
	res = addWithMask((BYTE *)imdWithMask + 2 * btsVect, (const BYTE *)product, (const BYTE *)rvectWithMask + 2 * btsVect, dimsB);

	res = bitAndWithMask((BYTE *)product, (const BYTE *)imdWithMask + 1 * btsVect, (const BYTE *)imdWithMask + 2 * btsVect, dimsB);
	res = addWithMask((BYTE *)imdWithMask + 3 * btsVect, (const BYTE *)product, (const BYTE *)rvectWithMask + 1 * btsVect, dimsB);

	/* Generate the correct order */
	/* d a b c */
	//3 0 1 2
	for (i = 0; i < MASKD; ++i){
		*(s4bRes + i * btsMat)= imdWithMask[3][i][0];
		*(s4bRes + i * btsMat + 1) = imdWithMask[3][i][1];

		*(s4bRes + i * btsMat + 2) = imdWithMask[0][i][0];
		*(s4bRes + i * btsMat + 3) = imdWithMask[0][i][1];

		*(s4bRes + i * btsMat + 4) = imdWithMask[1][i][0];
		*(s4bRes + i * btsMat + 5) = imdWithMask[1][i][1];

		*(s4bRes + i * btsMat + 6) = imdWithMask[2][i][0];
		*(s4bRes + i * btsMat + 7) = imdWithMask[2][i][1];

	}
	return res;
}

#else /* Unmask */
/* A 4-bit SBOX (unmask) */
static
Res sbox4b(
BYTE *s4bRes,
const BYTE *mats4b,
const int *dims4b// [4, DIM_L]
)
{
	if (dims4b[0] != 4) return RES_INVALID_DIMENSION;
	Res res = RES_OK;
	BYTE rvect[4][2] = { 0 };
	BYTE imd[4][2] = { 0 };
	//rvect = split(mat4b, 4, 1);
	int j;
	for (j = 0; j < 4; ++j){
		rvect[j][0] = mats4b[j * 2];
		rvect[j][1] = mats4b[j * 2 + 1];
	}

	// a = 1 x 2 + 0
	// b = 2 x a + 3
	// c = 1 x b + 2
	// d = b x c + 1


	BYTE product[2] = { 0 };
	int dimsB[4] = { 1, DIM_L, 1, DIM_L };
	int btsVect = 2;
	res = bitAnd(product, (const BYTE *)rvect + btsVect, (const BYTE *)rvect + 2 * btsVect, dimsB);
	res = add((BYTE *)imd, (const BYTE *)product, (const BYTE *)rvect, dimsB);

	res = bitAnd(product, (const BYTE *)rvect + btsVect * 2, (const BYTE *)imd, dimsB);
	res = add((BYTE *)imd + btsVect, (const BYTE *)product, (const BYTE *)rvect + btsVect * 3, dimsB);

	res = bitAnd(product, (const BYTE *)rvect + btsVect, (const BYTE *)imd + btsVect, dimsB);
	res = add((BYTE *)imd + 2 * btsVect, (const BYTE *)product, (const BYTE *)rvect + btsVect * 2, dimsB);

	res = bitAnd(product, (const BYTE *)imd[1], (const BYTE *)imd + btsVect * 2, dimsB);
	res = add((BYTE *)imd + 3 * btsVect, (const BYTE *)product, (const BYTE *)rvect + btsVect, dimsB);


	/* Generate the correct order */
	/* d a b c */ // 3012
	*(s4bRes) = imd[3][0];
	*(s4bRes + 1) = imd[3][1];

	*(s4bRes + 2) = imd[0][0];
	*(s4bRes + 3) = imd[0][1];

	*(s4bRes + 4) = imd[1][0];
	*(s4bRes + 5) = imd[1][1];

	*(s4bRes + 6) = imd[2][0];
	*(s4bRes + 7) = imd[2][1];

	return res;
}
#endif /* MASK */



/* LBOX */
#if MASK
#if DIVIDE
static
void lboxes(
BYTE **matsLin,
int direction /* direction == 0, means the left */
)
{
	BYTE *matTem = matsLin[0];
	matsLin[0] = multiply(matTem, TSlices[direction]);
	deMat(matTem);

	int indexOfMask;
	for (indexOfMask = 1; indexOfMask < MASKD; ++indexOfMask){
		matTem = matsLin[indexOfMask];
		matsLin[indexOfMask] = multiply(matTem, LSlices[direction]);
		deMat(matTem);
	}
}

#else /* DIVIDE, DIM_A == DIM_L || DIM_A == 0 */
static
Res  lboxes(
BYTE *matsLin
)
{
	Res res;
	int i;
	/* z[0] = T  *  x[0] */
	BYTE matTem[DIM_L] = { 0 };
	const int dims[4] = { DIM_S, DIM_L, DIM_L, DIM_L };

	memcpy(matTem, (const BYTE*)matsLin, DIM_L);
#if DIM_A
	res = multiply(matsLin, (const BYTE *)matTem, (const BYTE *)matT, dims);
#else /* DIM_A == 0 */
	res = multiply(matsLin, (const BYTE *)matTem, (const BYTE *)matL, dims);
#endif
	
	for (i = 1; i != MASKD; ++i){
		memcpy(matTem, (const BYTE*)matsLin + i * DIM_L, DIM_L);
		res = multiply(matsLin + i * DIM_L, (const BYTE *)matTem, (const BYTE *)matL, dims);
	}
	return res;
}
#endif /* DIVIDE */
#else /* Unmask  */
/* DIM_L-bit L-box (unmask) */
static
Res lboxes(
BYTE *lin
)
{
	Res res = RES_OK;
	const int dims[4] = { DIM_S, DIM_L, DIM_L, DIM_L };
	BYTE tem[DIM_L] = { 0 };
	memcpy(tem, lin, DIM_L);
	res = multiply(lin, tem, (const BYTE *)matL, dims);
	return res;
}
#endif /* MASK */

/* SBOX*/
#if MASK
/* DIM_S-bit S-box (masked)*/
static
Res  sboxes(
BYTE *matsSin,
const BYTE *keyRound
)
{
	int i;
	Res res = RES_OK;
	BYTE s4bRes[DIM_L][DIM_L] = { 0 };
	BYTE left[MASKD][DIM_L / 2] = { 0 };
	BYTE right[MASKD][DIM_L / 2] = { 0 };
	BYTE sum[MASKD][DIM_L / 2] = { 0 };
	
	BYTE *ptrOfL, *ptrOfR;
	int  btsHalf;
	int  theLast;
	const int dimsF[4] = { DIM_S / 2, DIM_L, DIM_S / 2, DIM_L };

	btsHalf = DIM_L / 2;
	theLast = MASKD - 1;

	for (i = 0; i < MASKD; ++i){
		memcpy(left + i * btsHalf, (const BYTE *)matsSin + i * btsHalf * 2, btsHalf);
		memcpy(right + i * btsHalf, (const BYTE *)matsSin + i *btsHalf * 2 + 1, btsHalf);
	}

	/* Combine a bigger sbox, from 4-bit to 8-bit */
	
	ptrOfL = (BYTE *)left, ptrOfR = (BYTE *)right;
	/* Feistel Struct Begins */
	for (i = 0; i < FEISTEL; ++i)
	{
		BYTE tem[DIM_L / 2] = { 0 };
		memcpy(tem, (const BYTE*)ptrOfL + theLast * btsHalf, btsHalf);
		/* add key_r to the last masked component */
		res = add(ptrOfL + theLast * btsHalf, (const BYTE *)tem, (const BYTE *)keyRound, dimsF);

		/* pass a 4-bit sbox */
		res = sbox4b((BYTE *)s4bRes, (const BYTE *)ptrOfL, dimsF);

		/* recover the left matrices */
		memcpy(ptrOfL + theLast * btsHalf, (const BYTE*)tem, btsHalf);

		/* do 'XOR' with the right matrix */
		res = addWithMask((BYTE *)sum, (const BYTE *)ptrOfR, (const BYTE *)s4bRes, dimsF);

		/* exchange  */
		ptrOfR = ptrOfL;
		ptrOfL = (BYTE *)sum;
	}

	/* Catenate those vectors to a matrix to return */
	for (i = 0; i < MASKD; ++i){
		memcpy(matsSin + i * btsHalf * 2, (const BYTE *)ptrOfL + i *btsHalf, btsHalf);
		memcpy(matsSin + i *btsHalf * 2 + 1, (const BYTE *)ptrOfR + i *btsHalf, btsHalf);
	}
	return res;
}

#else /* Unmask */

/* DIM_S-bit S-box (unmask)*/
static
Res sboxes(
BYTE *sin
)
{
	/* Split matrix to rowsUpper and rowsLower */

	BYTE s4bRes[DIM_L][DIM_L] = { 0 };
	BYTE leftT[DIM_L / 2] = { 0 };
	BYTE rightT[DIM_L / 2] = { 0 };
	BYTE sumT[DIM_L / 2] = { 0 };

	Res res = RES_OK;
	const int dimsF[4] = { DIM_S / 2, DIM_L, DIM_S / 2, DIM_L };
	int btsHalf = DIM_L / 2;

	memcpy(leftT , sin , btsHalf);
	memcpy(rightT, sin + btsHalf, btsHalf);
	BYTE *sum = sumT, *left = leftT, *right = rightT;

	/* Combine a bigger sbox, from 4-bit to 8-bit */
	int i;
	/* Feistel Struct */
	for (i = 0; i < FEISTEL; ++i)
	{
		/* Firstly, do 'XOR' with key_r */
		res = add(sum, left, keyR, dimsF);

		/* Secondly, through a 4-bit sbox */
		res = sbox4b(s4bRes, sum, dimsF);

		/* Then do 'XOR' with the right matrix */
		res = add(sum, right, s4bRes, dimsF);
	
		/* Finally, exchange each side */
		right = left;
		left = sum;
	}


	/* Catenate those vectors to a matrix to return */
	memcpy(sin, left , btsHalf);
	memcpy(sin + btsHalf * 2 + 1, right , btsHalf);

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
	res = multiply(matRight, (const BYTE *)matTransAs, (const BYTE *)matL, dims);
	res = multiply(matT, (const BYTE *)matInvAs, (const BYTE *)matRight, dims);
#else
	res = multiply(matRight, (const BYTE *)matTransA, (const BYTE *)matL, dims);
	res = multiply(matT, (const BYTE *)matInvA, (const BYTE *)matRight, dims);
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

#if DIVIDE
{
	/* Split to slices in vertical dimension */

	BYTE **keySlices = split(key, DIVIDE_PARTS, 2);
	BYTE **plainSlices = split(plain, DIVIDE_PARTS, 2);

	BYTE *cipherSlices[DIVIDE_PARTS];

	/* Encoded Plain */
	int indexOfSlices;
	int theLast = MASKD - 1;
	BYTE **maskedPlain[DIVIDE_PARTS] = { 0 };
	for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
		maskedPlain[indexOfSlices] = encode(plainSlices[indexOfSlices]);
		if (maskedPlain[indexOfSlices] == NULL) return NULL;

		/* Add Key */
		BYTE *matTem = maskedPlain[indexOfSlices][theLast];
		maskedPlain[indexOfSlices][theLast] = add(maskedPlain[indexOfSlices][theLast], keySlices[indexOfSlices]);
		deMat(matTem);
	}

	int indexOfRound;
	BYTE **matsMix = NULL;
	BYTE **matsTem = NULL;
	for (indexOfRound = 0; indexOfRound != ROUNDS; ++indexOfRound)
	{
		for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
			/* S-box */
			sboxes(maskedPlain[indexOfSlices], keyRoundSlices[indexOfSlices]);
			/* L-box */
			lboxes(maskedPlain[indexOfSlices], indexOfSlices);
			/* Mix the slices */
			if (!indexOfSlices) matsMix = maskedPlain[indexOfSlices];
			else {
				matsTem = matsMix;
				matsMix = addWithMask(maskedPlain[indexOfSlices], matsTem);
				deMats(maskedPlain[indexOfSlices], MASKD);
				deMats(matsTem, MASKD);
			}
			//maskedPlain[indexOfSlices] = (BYTE **)malloc(MASKD * sizeof(BYTE *));
		}  
		/* Split it */
		int indexOfMask;
		for (indexOfMask = 0; indexOfMask != MASKD; ++indexOfMask){
			matsTem = split(matsMix[indexOfMask], DIVIDE_PARTS, 2);
			for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
				maskedPlain[indexOfSlices][indexOfMask] = matsTem[indexOfSlices];
			}			
		}


		/* Add Key And Round Constant */
		for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
			BYTE *roundKeySlice = add(rdConstSlices[indexOfRound][indexOfSlices], keySlices[indexOfSlices]);
			BYTE *matTem = maskedPlain[indexOfSlices][theLast];
			maskedPlain[indexOfSlices][theLast] = add(roundKeySlice, matTem);
			deMat(matTem);
			deMat(roundKeySlice);
		}

	}

	for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
		cipherSlices[indexOfSlices] = decode(maskedPlain[indexOfSlices]);
	}

	BYTE *cipher = cat(cipherSlices, DIVIDE_PARTS, 2);

	// dePostCal();
	return cipher;
}

#else /* DIVIDE */
{
	Res res = RES_OK;
	BYTE matMasked[MASKD][DIM_L] = { 0 };
	BYTE tem[DIM_L] = { 0 };
	int theLast;
	const int dimsA[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	int indexOfRound;
	int btsMat;
	if (plain == NULL || key == NULL) return RES_INVALID_POINTER;

	/* Encoded Plain */

	res = encode((BYTE *)matMasked, plain);

	/* Add Key */
	theLast = MASKD - 1;
	btsMat = DIM_S * bytesOfRow(DIM_L);
	

	memcpy(tem, (const BYTE *)(matMasked + theLast * btsMat), DIM_L);
	res = add((BYTE *)(matMasked + theLast * btsMat), key, (const BYTE *)tem, dimsA);	
	for (indexOfRound = 0; indexOfRound < ROUNDS; ++indexOfRound)
	{
		BYTE roundK[DIM_L] = { 0 };

		res = sboxes((BYTE *)matMasked, (BYTE *)keyR);
		res = lboxes((BYTE *)matMasked);
		res = add((BYTE *)roundK, (const BYTE *)rdConst, key, dimsA);

		/* Add Key And Round Constant */
		memcpy(tem, (const BYTE *)(matMasked + theLast * btsMat), DIM_L);
		res = add((BYTE *)(matMasked + theLast * btsMat), (const BYTE *)tem, (const BYTE *)roundK, dimsA);

	}
	/* Decode Cipher */
	res = decode(cipher, (const BYTE *)matMasked);
	return res;

}
#endif /* DIVIDE_PARTS */

#else /* Unmask */
	/* Encryption begins */
{
	if (plain == NULL || key == NULL || ROUNDS < 0) return RES_INVALID_POINTER;

	BYTE  tem[DIM_L] = { 0 };
	Res res = RES_OK;
	int btsMat = DIM_S * bytesOfRow(DIM_L);
	const int dimsA[4] = { DIM_S, DIM_L, DIM_S, DIM_L };
	res = add(cipher, plain, key, dimsA);
	
	int round_i;
	for (round_i = 0; round_i != ROUNDS; ++round_i)
	{
		res = sboxes(cipher);
		res = lboxes(cipher);
		memcpy(tem, cipher, DIM_L);
		res = add(tem,cipher, key, dimsA);

		res = add(cipher, tem, rdConst, dimsA);
	}

	return res;
}
#endif /* MASK */

void encrypto_fixed(){
	Res res = RES_OK;
	BYTE cipher[DIM_L] = { 0 };

	const BYTE plainT[DIM_L] = { 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf, 0xa0 };
	const BYTE cipherK[DIM_L] = { 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf, 0xb0 };

	res = encrypto(cipher, plainT, cipherK);
}