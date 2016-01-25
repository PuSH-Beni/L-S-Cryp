//
//  LSbox.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#include "LSbox.h"

/*=========================================================*/
/*   MARK: Declarition     */
/*=========================================================*/

/* round_constant , key_round, matrix L*/
static
Mat  *rdConst[ROUNDS],  *key_r, *matL;

#if MASK

#if DIM_A == 4
static
Mat  *rdCLeft[ROUNDS], *rdCRight[ROUNDS], *Lleft, *Lright, *Tleft, *Tright, **keyRoundLR;
#endif

#if DIM_A != 0
extern
Mat *matA, *matInvA, *matTransA;
static
Mat *matT;
#endif

#endif /* MASK */

BYTE matLV[] = MAT_LV;
BYTE keyRV[] = KEY_RV;
BYTE rdConst1V[] = CONSTR1;
BYTE rdConst2V[] = CONSTR2;
BYTE rdConst3V[] = CONSTR3;

/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/

#if MASK
#if DIM_A == 4
static
void splitHorizonParts(
            const Mat **matsO,
            Mat **matLeft,
            Mat **matRight,
            int parts
            )
{
    int index;
    Mat ***matsRet = malloc(parts * sizeof(Mat **));
    for (index = 0; index < parts; ++index){
        matsRet[index] = split(matsO[index], 2, 2);
        matLeft[index] = matsRet[index][0];
        matRight[index] = matsRet[index][1];
    }
}
#endif
#endif /* MASK */



#if MASK
/* A 4-bit sbox (masked)*/
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

    product = bitAndWithMask(rvectWithMask[1], rvectWithMask[2]);
    imdWithMask[0] = addWithMask(product, rvectWithMask[0]);
    deMats(product);

    product = bitAndWithMask(rvectWithMask[2], imdWithMask[0]);
    imdWithMask[1] = addWithMask(product, rvectWithMask[3]);
    deMats(product);

    product = bitAndWithMask(rvectWithMask[1], imdWithMask[1]);
    imdWithMask[2] = addWithMask(product, rvectWithMask[2]);
    deMats(product);

    product = bitAndWithMask(imdWithMask[1], imdWithMask[2]);
    imdWithMask[3] = addWithMask(product, rvectWithMask[1]);
    deMats(product);

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
        retMat[indexOfMask] = cat(matsTem[indexOfMask], 4, 1);
    }


    /* Deallocate all  */
    int i;
    for(i = 0; i < 4; ++i)
    {
        deMats(imdWithMask[i]);
    }

    return retMat;
}

#else /* Unmask */
/* A 4-bit sbox (unmask) */
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
    
    product = bitAnd(rvect[1], rvect[2]);
    imd[0] = add(product, rvect[0]);
    deMat(product);
    
    product = bitAnd(rvect[2], imd[0]);
    imd[1] = add(product, rvect[3]);
    deMat(product);
    
    product = bitAnd(rvect[1], imd[1]);
    imd[2] = add(product, rvect[2]);
    deMat(product);
    
    product = bitAnd(imd[1], imd[2]);
    imd[3] = add(product, rvect[1]);
    deMat(product);


    /* Generate the correct order */
    /* d a b c */
    Mat *ordered[] = {imd[3], imd[0], imd[1], imd[2]};

    Mat *retMat;
    retMat =  cat(ordered, 4, 1);

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


/* Using matT, matL */
/* A Lboxes (masked) */
#if MASK
#if DIM_A == 4
static
void lboxes(
            Mat **matsLin,
            int direction /* direction == 0, means the left */
            )
{
    Mat *matTem = matsLin[0];
    int i;
	if (!direction) {
        matsLin[0] = multiply(matTem,Tleft);
    }else{
        matsLin[0] = multiply(matTem,Tright);
        /*for (i = 0; i < matsLin[0]->dim_row; ++i) {
            matsLin[0]->vect[i] = matsLin[0]->vect[i] >> 4;
        }*/
    }
    deMat(matTem);
    matsLin[0]->dim_col = 8;


    for(i = 1; i < MASKD; ++i){
        matTem = matsLin[i];
        if (!direction) {
            matsLin[i] = multiply(matTem,Lleft);
        }else{
            matsLin[i] = multiply(matTem,Lright);
        }
        deMat(matTem);
    }
}

#else /* DIM_A ==8 or 0 */
static
void lboxes(
            Mat **matsLin
            )
{
    /* z_1 = T  x x_1 */
    Mat *matTem = matsLin[0];
#if DIM_A == 8
    matsLin[0] = multiply(matTem,matT);
#else /* DIM_A == 0 */
    matsLin[0] = multiply(matTem, matL);
#endif
    deMat(matTem);

    int i;
    for(i = 1; i < MASKD; ++i){
        matTem = matsLin[i];
        matsLin[i] = multiply(matTem,matL);
        deMat(matTem);
    }
}
#endif /* DIM_A */

#else /* Unmask  */
/* DIM_L-bit L-box (unmask) */
static
Mat *lboxes(
            const Mat *lin
            )
{
    Mat *retMat;
    retMat = multiply(lin,matL);
    return retMat;
}
#endif /* MASK */


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
        rows[indexOfMask] = split(matsSin[indexOfMask], 2, 1);
        
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
        ptrOfL[theLast] = add(ptrOfL[theLast], keyRound);
        
        
        /* Secondly, through a 4-bit sbox */
        fout4b = sbox4b(ptrOfL);
        
        /* Then do 'XOR' with the right matrix */
        /* Get the NEXT LEFT part  */
        
        sum = addWithMask(ptrOfR, fout4b);
        
        deMats(fout4b);
        deMats(ptrOfR);
        
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
        matsSin[indexOfMask] = cat(rows[indexOfMask], 2, 1);

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
static
void newPreCal()

{

    key_r = newMat(DIM_S / 2, DIM_L, keyRV, 0x03);
    rdConst[0] = newMat(DIM_S, DIM_L, rdConst1V, 0x03);
    rdConst[1] = newMat(DIM_S, DIM_L, rdConst2V, 0x03);
    rdConst[2] = newMat(DIM_S, DIM_L, rdConst3V, 0x03);
    matL = newMat(DIM_S, DIM_L, matLV, 0x03);


#if MASK

#if DIM_A == 4
    /* Get Matrix T-LR and L-LR*/
    Mat **LLR = split(matL, 2, 2);
    Lleft = LLR[0];
    Lright = LLR[1];

    /* Get Matrix T  */ 
    //Mat *matDoubleInv = cat(matsTem, 2, 2);
	/*  matDoubleInv looks like: 
    	| inv(A), 0      |
	    |      0, inv(A) | */
	Mat *matDoubleInv = newMat(DIM_A*2, DIM_A*2, NULL, 0x00);
	
	memmove(matDoubleInv->vect, matInvA->vect, DIM_A * sizeof(BYTE));
	memmove(matDoubleInv->vect + DIM_A, matInvA->vect, DIM_A * sizeof(BYTE));
	/* Shift the bits of the right-lower matrix */
	int i;
	for (i = DIM_A; i < matDoubleInv->dim_row; ++i){
		matDoubleInv->vect[i] = matDoubleInv->vect[i] >> 4;
	}

    Mat *matRightPart = multiply(matTransA, Lleft);
    Tleft= multiply(matDoubleInv, matRightPart);
    deMat(matRightPart);

    matRightPart = multiply(matTransA, Lright);
    Tright = multiply(matDoubleInv, matRightPart);

    deMat(matRightPart);
    deMat(matDoubleInv);

    splitHorizonParts(rdConst, rdCLeft, rdCRight, ROUNDS);
    keyRoundLR = split(key_r, 2, 2);

#elif DIM_A == 8
    /* Get Matrix T  */
    Mat *matRight = multiply(matTransA, matL);
    matT = multiply(matInvA, matRight);

    deMat(matRight);

#endif /* DIM_A */

#endif /* MASK */
}



/* After encryption, deconstruct those matrices */
static
void dePostCal()
{

    deMat(key_r);
    int i;
    for (i = 0; i < ROUNDS; ++i){
        deMat(rdConst[i]);
    }
#if MASK 
#if DIM_A /* USING A */
    deMat(matT);
#endif
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
#if DIM_A == 4
{
    setup();
    newPreCal();
    /* Split to 2 parts in vertical dimension */
    //Mat *plainLeft[MASKD], *plainRight[MASKD];
    Mat **keyLR = split(key, 2, 2);
    Mat **plainLR = split(plain, 2, 2);

    Mat *ciphers[2];

    /* Encoded Plain */
    Mat **maskedL = encode(plainLR[0]);
    Mat **maskedR = encode(plainLR[1]);
    if (maskedL == NULL) return NULL;
    if (maskedR == NULL) return NULL;

    /* Add Key */
    int theLast = MASKD - 1;
    Mat *matTem = maskedL[theLast];
    maskedL[theLast] = add(maskedL[theLast], keyLR[0]);
    deMat(matTem);

    matTem = maskedR[theLast];
    maskedR[theLast] = add(maskedR[theLast], keyLR[1]);
    deMat(matTem);

    int indexOfRound;
    for (indexOfRound = 0; indexOfRound < ROUNDS; ++indexOfRound)
    {
        /* S-box */
        sboxes(maskedL, keyRoundLR[0]);
        sboxes(maskedR, keyRoundLR[1]);

        /* L-box */
        lboxes(maskedL, 0);
        lboxes(maskedR, 1);

        /* Mix the left and the right */
        Mat **matsMix = addWithMask(maskedL, maskedR);

        deMats(maskedL);
        deMats(maskedR);
        
        splitHorizonParts(matsMix, maskedL, maskedR, MASKD);
        
        //Res res =refreshing(maskedL);
        //res =refreshing(maskedR);

        Mat *roundKeyL = add(rdCLeft[indexOfRound], keyLR[0]);
        Mat *roundKeyR = add(rdCRight[indexOfRound], keyLR[1]);

        /* Add Key And Round Constant */
        matTem = maskedL[theLast];
        maskedL[theLast] = add(maskedL[theLast], roundKeyL);

        deMat(matTem);
        deMat(roundKeyL);

        matTem = maskedR[theLast];
        maskedR[theLast] = add(maskedR[theLast], roundKeyR);

        deMat(matTem);
        deMat(roundKeyR);
    }

    ciphers[0] = decode(maskedL);
    ciphers[1] = decode(maskedR);

    Mat *cipher = cat(ciphers, 2, 2);

    dePostCal();
    return cipher;
}

#else /* DIM_A == 8 or 0 */
{
    if (plain == NULL || key == NULL || ROUNDS < 0) return NULL;
    Mat  *matRoundKey,  *matTem;

#if DIM_A /* DIM_A != 0 */
    setup();
#endif

    newPreCal();

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
        
        //Res res =refreshing(matsMasked);

        matRoundKey = add(rdConst[indexOfRound], key);

        /* Add Key And Round Constant */
        matTem = matsMasked[theLast];
        matsMasked[theLast] = add(matsMasked[theLast], matRoundKey);

        deMat(matTem);
        deMat(matRoundKey);

    }

    /* Decode Cipher */
    Mat *cipher = decode(matsMasked);

    dePostCal();

    return cipher;

}
#endif /* DIM_A  */

#else /* Unmask */
/* Encryption begins */
{
    if(plain == NULL || key == NULL || ROUNDS < 0) return NULL;

    Mat  *roundIn, *sum, *sout, *lout;	
    newPreCal();

    roundIn = add(plain, key);
    int round_i;
    for(round_i = 0; round_i < ROUNDS; ++round_i)
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

    dePostCal();
    return roundIn;
}
#endif /* MASK */
