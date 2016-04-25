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

/* round_constant , key_round, matrix L*/
static 
Mat  *rdConst[ROUNDS],  *key_r, *matL;

#if MASK
#if DIVIDE
static 
Mat **keyRoundSlices, **rdConstSlices[ROUNDS], **LSlices, *TSlices[DIVIDE_PARTS];
#endif

#if DIM_A
extern
Mat *matA, *matInvA, *matTransA;
#if !DIVIDE && DIM_L == 16
extern
Mat *matAs, *matInvAs, *matTransAs;
#endif

static
Mat  *matT;
#endif /* DIM_A */

#endif /* MASK */
/********************************************
 *  static constants
********************************************/
BASE matLV[L_SIZE] = MAT_LV;
BASE keyRV[KEY_SIZE] = KEY_RV;
BASE rdConst1V[CONST_SIZE] = CONSTR1;
BASE rdConst2V[CONST_SIZE] = CONSTR2;
BASE rdConst3V[CONST_SIZE] = CONSTR3;
BASE rdConst4V[CONST_SIZE] = CONSTR4;
BASE rdConst5V[CONST_SIZE] = CONSTR5;
BASE rdConst6V[CONST_SIZE] = CONSTR6;
BASE rdConst7V[CONST_SIZE] = CONSTR7;
BASE rdConst8V[CONST_SIZE] = CONSTR8;

BASE rdConst9V[CONST_SIZE] = CONSTR9;
BASE rdConst10V[CONST_SIZE] = CONSTR10;
BASE rdConst11V[CONST_SIZE] = CONSTR11;
BASE rdConst12V[CONST_SIZE] = CONSTR12;
#if ROUNDS > 12
BASE rdConst13V[CONST_SIZE] = CONSTR13;
BASE rdConst14V[CONST_SIZE] = CONSTR14;
BASE rdConst15V[CONST_SIZE] = CONSTR15;
BASE rdConst16V[CONST_SIZE] = CONSTR16;
#endif


/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/

/* 4-bit SBOX */
#if MASK
static
Mat **sbox4b(
             const Mat **mats4b
             )
{
    if (mats4b[0]->dim_row != 4) return NULL;

    Mat **matsTem[MASKD];
    Mat *rvectWithMask[4][MASKD], **imdWithMask[4];

    int indexOfMask;
    for (indexOfMask = 0; indexOfMask < MASKD; ++indexOfMask){
        matsTem[indexOfMask] = split(mats4b[indexOfMask], 4, 1);
        int i;
        for (i = 0; i < 4; ++i){
            rvectWithMask[i][indexOfMask] = matsTem[indexOfMask][i];
        }
    }

    // a = 1 x 2 + 0
    // b = 2 x a + 3
    // c = 1 x b + 2
    // d = b x c + 1
    Mat **product;

    product = bitAndWithMask((const Mat **)rvectWithMask[1], (const Mat **)rvectWithMask[2]);
    imdWithMask[0] = addWithMask((const Mat **)product, (const Mat **)rvectWithMask[0]);
    deMats(product, MASKD);

    product = bitAndWithMask((const Mat **)rvectWithMask[2], (const Mat **)imdWithMask[0]);
    imdWithMask[1] = addWithMask((const Mat **)product, (const Mat **)rvectWithMask[3]);
    deMats(product, MASKD);

    product = bitAndWithMask((const Mat **)rvectWithMask[1], (const Mat **)imdWithMask[1]);
    imdWithMask[2] = addWithMask((const Mat **)product, (const Mat **)rvectWithMask[2]);
    deMats(product, MASKD);

    product = bitAndWithMask((const Mat **)imdWithMask[1], (const Mat **)imdWithMask[2]);
    imdWithMask[3] = addWithMask((const Mat **)product, (const Mat **)rvectWithMask[1]);
    deMats(product, MASKD);

    /* Generate the correct order */
    /* d a b c */
    Mat **ordered[] = { imdWithMask[3], imdWithMask[0], imdWithMask[1], imdWithMask[2] };

    Mat **retMat = (Mat **)malloc(MASKD * sizeof(Mat *));
    for (indexOfMask = 0; indexOfMask < MASKD; ++indexOfMask){
        int i;
        for (i = 0; i < 4; ++i){
            deMat(matsTem[indexOfMask][i]);
            matsTem[indexOfMask][i] = ordered[i][indexOfMask];
        }
        retMat[indexOfMask] = cat((const Mat **)matsTem[indexOfMask], 4, 1);
    }


    /* Deallocate all  */
    int i;
    for(i = 0; i < 4; ++i)
    {
        deMats(imdWithMask[i], MASKD);
    }

    return retMat;
}

#else /* Unmask */
/* A 4-bit SBOX (unmask) */
static
Mat *sbox4b(
            const Mat *mat4b
            )
{
    if(mat4b->dim_row != 4) return NULL;
    
    Mat **rvect;
    rvect = split(mat4b, 4, 1);
    
    // a = 1 x 2 + 0
    // b = 2 x a + 3
    // c = 1 x b + 2
    // d = b x c + 1
    
    Mat *imd[4];
    Mat *product;
    
    product = bitAnd((const Mat *)rvect[1], (const Mat *)rvect[2]);
    imd[0] = add((const Mat *)product, (const Mat *)rvect[0]);
    deMat(product);
    
    product = bitAnd((const Mat *)rvect[2], (const Mat *)imd[0]);
    imd[1] = add((const Mat *)product, (const Mat *)rvect[3]);
    deMat(product);
    
    product = bitAnd((const Mat *)rvect[1], (const Mat *)imd[1]);
    imd[2] = add((const Mat *)product, (const Mat *)rvect[2]);
    deMat(product);
    
    product = bitAnd((const Mat *)imd[1], (const Mat *)imd[2]);
    imd[3] = add((const Mat *)product, (const Mat *)rvect[1]);
    deMat(product);


    /* Generate the correct order */
    /* d a b c */
    Mat *ordered[] = {imd[3], imd[0], imd[1], imd[2]};

    Mat *retMat;
    retMat =  cat((const Mat**)ordered, 4, 1);

    /* Deallocate all  */
    int i;
    for(i = 0; i < 4; ++i)
    {
        deMat(imd[i]);
        deMat(rvect[i]);
    }

    return retMat;
}
#endif /* MASK */



/* LBOX */
#if MASK
#if DIVIDE
static
void lboxes(
            Mat **matsLin,
            int direction /* direction == 0, means the left */
            )
{
    Mat *matTem = matsLin[0];
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
void lboxes(
            Mat **matsLin
            )
{
    /* z[0] = T  *  x[0] */
    Mat *matTem = matsLin[0];
#if DIM_A
    matsLin[0] = multiply((const Mat *)matTem, (const Mat *)matT);
#else /* DIM_A == 0 */
    matsLin[0] = multiply((const Mat *)matTem, (const Mat *)matL);
#endif
    deMat(matTem);

    int i;
    for(i = 1; i != MASKD; ++i){
        matTem = matsLin[i];
        matsLin[i] = multiply((const Mat *)matTem, (const Mat *)matL);
        deMat(matTem);
    }
}
#endif /* DIVIDE */
#else /* Unmask  */
/* DIM_L-bit L-box (unmask) */
static
Mat *lboxes(
            const Mat *lin
            )
{
    Mat *retMat;
    retMat = multiply(lin, (const Mat *)matL);
    return retMat;
}
#endif /* MASK */

/* SBOX*/
#if MASK
/* DIM_S-bit S-box (masked)*/
static
void sboxes(
            Mat **matsSin,
			Mat *keyRound
            )
{
    /* Split matrices to rowsUpper and rowsLower */
    Mat **rows[MASKD], *left[MASKD], *right[MASKD];
    
    int indexOfMask;
    for (indexOfMask = 0; indexOfMask < MASKD; ++indexOfMask){
        rows[indexOfMask] = split((const Mat *)matsSin[indexOfMask], 2, 1);
        
        /* Get the Left  and the Right parts */
        left[indexOfMask] = rows[indexOfMask][0];
        right[indexOfMask] = rows[indexOfMask][1];
    }
    
    
    Mat **ptrOfL = left, **ptrOfR = right;
    
    /* Combine a bigger sbox, from 4-bit to 8-bit */
    int i, theLast = MASKD - 1;
    Mat **sum = NULL, **fout4b = NULL, *matOlder = NULL;

    /* Feistel Struct Begins */
    for (i = 0; i < FEISTEL; ++i)
    {
        /* Firstly, add key_r to the last masked component */

        matOlder = ptrOfL[theLast];
        ptrOfL[theLast] = add((const Mat *)ptrOfL[theLast], (const Mat *)keyRound);
        
        
        /* Secondly, through a 4-bit sbox */
        fout4b = sbox4b((const Mat **)ptrOfL);
        
        /* Then do 'XOR' with the right matrix */
        /* Get the NEXT LEFT part  */
        
        sum = addWithMask((const Mat **)ptrOfR, (const Mat **)fout4b);
        
        deMats(fout4b, MASKD);
        deMats(ptrOfR, MASKD);
        
        /* Finally, recover the left matrices */
        deMat(ptrOfL[theLast]);
        ptrOfL[theLast] = matOlder;
        ptrOfR = ptrOfL;
        ptrOfL = sum;
    }
    
    
    /* Catenate those vectors to a matrix to return */
        
    for (indexOfMask = 0; indexOfMask < MASKD; ++indexOfMask){
        rows[indexOfMask][0] = ptrOfL[indexOfMask];
        rows[indexOfMask][1] = ptrOfR[indexOfMask];
        
        deMat(matsSin[indexOfMask]);
        matsSin[indexOfMask] = cat((const Mat **)rows[indexOfMask], 2, 1);

        /* Deallocate all  */
        deMat(ptrOfL[indexOfMask]);
        deMat(ptrOfR[indexOfMask]);
        free(rows[indexOfMask]);
    }

}

#else /* Unmask */

/* DIM_S-bit S-box (unmask)*/
static
Mat *sboxes(
            const Mat *sin
            )
{
    /* Split matrix to rowsUpper and rowsLower */
    Mat **rows;
    rows = split(sin, 2, 1);


    /* Get the Left one and Right one */
    Mat *left, *right;
    left = rows[0];
    right = rows[1];

    /* Combine a bigger sbox, from 4-bit to 8-bit */
    int i;
    Mat *sum, *fout4b;
    /* Feistel Struct */
    for(i = 0; i < FEISTEL; ++i)
    {
        /* Firstly, do 'XOR' with key_r */
        sum = add(left, key_r);

        /* Secondly, through a 4-bit sbox */
        fout4b = sbox4b(sum);
        deMat(sum);

        /* Then do 'XOR' with the right matrix */
        sum = add(right, fout4b);
        deMat(fout4b);
        deMat(right);

        /* Finally, exchange each side */
        right = left;
        left = sum;
    }


    /* Catenate those vectors to a matrix to return */
    rows[0] = left;
    rows[1] = right;
    Mat *retMat;
    retMat = cat(rows, 2, 1);

    /* Deallocate all  */
    deMat(left);
    deMat(right);
    free(rows);

    return retMat;
}
#endif /* MASK */



/* Before encryption, do some pre-work to get the constant matrices */
void newPreCal()

{
	

	key_r = newMat((DIM_S / 2), DIM_L, keyRV, 0x03);
	matL = newMat(DIM_L, DIM_L, matLV, 0x03);

	rdConst[0] = newMat(DIM_S, DIM_L, rdConst1V, 0x03);
	rdConst[1] = newMat(DIM_S, DIM_L, rdConst2V, 0x03);
	rdConst[2] = newMat(DIM_S, DIM_L, rdConst3V, 0x03);
	rdConst[3] = newMat(DIM_S, DIM_L, rdConst4V, 0x03);
	rdConst[4] = newMat(DIM_S, DIM_L, rdConst5V, 0x03);
	rdConst[5] = newMat(DIM_S, DIM_L, rdConst6V, 0x03);
	rdConst[6] = newMat(DIM_S, DIM_L, rdConst7V, 0x03);
	rdConst[7] = newMat(DIM_S, DIM_L, rdConst8V, 0x03);

	rdConst[8] = newMat(DIM_S, DIM_L, rdConst9V, 0x03);
	rdConst[9] = newMat(DIM_S, DIM_L, rdConst10V, 0x03);
	rdConst[10] = newMat(DIM_S, DIM_L, rdConst11V, 0x03);
	rdConst[11] = newMat(DIM_S, DIM_L, rdConst12V, 0x03);
#if ROUNDS > 12
	rdConst[12] = newMat(DIM_S, DIM_L, rdConst13V, 0x03);
	rdConst[13] = newMat(DIM_S, DIM_L, rdConst14V, 0x03);
	rdConst[14] = newMat(DIM_S, DIM_L, rdConst15V, 0x03);
	rdConst[15] = newMat(DIM_S, DIM_L, rdConst16V, 0x03);
#endif

#if MASK

#if DIVIDE
	LSlices = split(matL, DIVIDE_PARTS, 2);
	/* Get Matrix T  */
	/*  matDoubleInv looks like:
		| inv(A), 0      |
		|      0, inv(A) | */

	Mat *matInvAs = newMat(DIM_A * DIVIDE_PARTS, DIM_A * DIVIDE_PARTS, NULL, 0x00);

	int btsOfRow = bytesOfRow(matInvAs->dim_col);
	int btsOfMat = matInvA->dim_row * bytesOfRow(matInvA->dim_col);
	int indexOfSlices;
#if DIM_A == 4  && DIM_L == 8	
	/* particular case*/
	memmove(matInvAs->vect, matInvA->vect, btsOfMat * sizeof(BASE));
	memmove(matInvAs->vect + DIM_A, matInvA->vect, btsOfMat * sizeof(BASE));
	/* Shift  bits  */
	int i;
	for(i = DIM_A; i != DIM_L; ++i) matInvAs->vect[i] >>= 4;				
#else 
	/* general cases */
	int indexOfByte, indexOfDest = 0;
	for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
#if DIM_A == 4
		BASE oddFlag = (BASE)indexOfSlices & 0x01; 
#endif
		for (indexOfByte = 0; indexOfByte != btsOfMat; ++indexOfByte){
#if DIM_A == 4
			matInvAs->vect[indexOfDest] = oddFlag ? matInvA->vect[indexOfByte] >> 4
													 : matInvA->vect[indexOfByte];
#else
			matInvAs->vect[indexOfDest] = matInvA->vect[indexOfByte];
#endif
			indexOfDest += btsOfRow;
		}
#if DIM_A == 4
		indexOfDest += oddFlag ? 1  : 0;
#else
		indexOfDest += 1;
#endif /* DIM_A == 4*/
	}
#endif/* get matInvAs */
		
	for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
		Mat *matRightPart = multiply(matTransA, LSlices[indexOfSlices]);
		TSlices[indexOfSlices] = multiply(matInvAs, matRightPart);
		deMat(matRightPart);
	}
	deMat(matInvAs);

	int indexOfRounds;
	for (indexOfRounds = 0; indexOfRounds != ROUNDS; ++indexOfRounds)
	{
		rdConstSlices[indexOfRounds] = split(rdConst[indexOfRounds], DIVIDE_PARTS, 2);
	}
	

	
#elif DIM_A /* !DIVIDE && DIM_A */
    /* Get Matrix T  */
#if DIM_L == 16
    Mat *matRight = multiply((const Mat *)matTransAs, (const Mat *)matL);
    matT = multiply((const Mat *)matInvAs, (const Mat *)matRight);
	//keyRoundSlices = split(key_r, DIVIDE_PARTS, 2);
    deMat(matRight);
#else
	Mat *matRight = multiply((const Mat *)matTransA, (const Mat *)matL);
	matT = multiply((const Mat *)matInvA, (const Mat *)matRight);

	deMat(matRight);
#endif

#endif /* DIVIDE */

#endif /* MASK */
}



/* After encryption, deconstruct those matrices */
void dePostCal()
{

    deMat(( Mat *)key_r);
	deMats((Mat **)rdConst, ROUNDS);
#if MASK  &&  DIM_A
/* USING A */
    deMat((Mat *)matT);
#endif

}



/*=========================================================*/
/*   MARK: Public   Functions      */
/*=========================================================*/

Mat *encrypto(
                const Mat *plain,
                const Mat *key
            )
#if MASK

#if DIVIDE
{
    /* Split to slices in vertical dimension */

    Mat **keySlices = split(key, DIVIDE_PARTS, 2);
    Mat **plainSlices = split(plain, DIVIDE_PARTS, 2);

    Mat *cipherSlices[DIVIDE_PARTS];

    /* Encoded Plain */
	int indexOfSlices;
	int theLast = MASKD - 1;
	Mat **maskedPlain[DIVIDE_PARTS] = { 0 };
	for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
		maskedPlain[indexOfSlices] = encode(plainSlices[indexOfSlices]);
		if (maskedPlain[indexOfSlices] == NULL) return NULL;
		
		/* Add Key */
		Mat *matTem = maskedPlain[indexOfSlices][theLast];
		maskedPlain[indexOfSlices][theLast] = add(maskedPlain[indexOfSlices][theLast], keySlices[indexOfSlices]);
		deMat(matTem);
	}
    
    int indexOfRound;
	Mat **matsMix = NULL;
	Mat **matsTem = NULL;
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
			//maskedPlain[indexOfSlices] = (Mat **)malloc(MASKD * sizeof(Mat *));
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
			Mat *roundKeySlice = add(rdConstSlices[indexOfRound][indexOfSlices], keySlices[indexOfSlices]);
			Mat *matTem = maskedPlain[indexOfSlices][theLast];
			maskedPlain[indexOfSlices][theLast] = add(roundKeySlice, matTem);
			deMat(matTem);
			deMat(roundKeySlice);
		}

    }

	for (indexOfSlices = 0; indexOfSlices != DIVIDE_PARTS; ++indexOfSlices){
		cipherSlices[indexOfSlices] = decode(maskedPlain[indexOfSlices]);
	}

	Mat *cipher = cat(cipherSlices, DIVIDE_PARTS, 2);

   // dePostCal();
    return cipher;
}

#else /* DIVIDE */
{
    if (plain == NULL || key == NULL || ROUNDS < 0) return NULL;
    Mat  *matRoundKey,  *matTem;

    /* Encoded Plain */
    Mat **matsMasked = encode(plain);
    if (matsMasked == NULL) return NULL;
    /* Add Key */
    int theLast = MASKD - 1;
    matTem = matsMasked[theLast];
    matsMasked[theLast] = add(matsMasked[theLast], key);
    deMat(matTem);

    int indexOfRound;
	for (indexOfRound = 0; indexOfRound < ROUNDS; ++indexOfRound)
	{
		sboxes(matsMasked, key_r);
        lboxes(matsMasked);      
		matRoundKey = add((const Mat *)rdConst[indexOfRound], (const Mat *)key);
       
		/* Add Key And Round Constant */
        matTem = matsMasked[theLast];
        matsMasked[theLast] = add((const Mat *)matsMasked[theLast], (const Mat *)matRoundKey);

        deMat(matTem);
        deMat(matRoundKey);

    }

    /* Decode Cipher */
    Mat *cipher = decode(matsMasked);

    return cipher;

}
#endif /* DIVIDE_PARTS */

#else /* Unmask */
/* Encryption begins */
{
    if(plain == NULL || key == NULL || ROUNDS < 0) return NULL;

    Mat  *roundIn, *sum, *sout, *lout;	

    roundIn = add(plain, key);
    int round_i;
    for(round_i = 0; round_i != ROUNDS; ++round_i)
    {
        sout = sboxes(roundIn);
        lout = lboxes(sout);
        deMat(roundIn);
        deMat(sout);

        sum = add(lout, key);
        deMat(lout);

        roundIn = add(sum, rdConst[round_i]);
        deMat(sum);
    }

    return roundIn;
}
#endif /* MASK */
