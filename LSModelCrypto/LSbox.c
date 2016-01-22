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
static
Mat  *rdConst[ROUNDS],  *key_r, *matL;
#if MASK && A_USING
static
Mat *matT;
#endif
BYTE matLV[] = MAT_LV;
BYTE keyRV[] = KEY_RV;
BYTE rdConst1V[] = CONSTR1;
BYTE rdConst2V[] = CONSTR2;
BYTE rdConst3V[] = CONSTR3;

/*=========================================================*/
/*   MARK: Private   Functions      */
/*=========================================================*/

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
#else
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
#endif

#if MASK
/* Using matT, matL */
/* A Lboxes (masked) */
static
void lboxes(
        Mat **matsLin
        )
{
	/* z_1 = T  x x_1 */
	Mat *matTem = matsLin[0];
#if A_USING
    matsLin[0] = multiply(matTem,matT);	
#else 
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

#else
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
#endif


#if MASK
/* DIM_S-bit S-box (masked)*/
static
void sboxes(
	Mat **matsSin
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
		ptrOfL[theLast] = add(ptrOfL[theLast], key_r);
		

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
#else

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
#endif


/* Before encryption, do some pre-work to get the constant matrices */
static
void newPreCal()

{

	key_r = newMat(DIM_S / 2, DIM_L, keyRV, 0x03);
	rdConst[0] = newMat(DIM_S, DIM_L, rdConst1V, 0x03);
	rdConst[1] = newMat(DIM_S, DIM_L, rdConst2V, 0x03);
	rdConst[2] = newMat(DIM_S, DIM_L, rdConst3V, 0x03);
	matL = newMat(DIM_S, DIM_L, matLV, 0x03);

#if MASK && A_USING   
	/* Get Matrix T  */
	Mat *matRight = multiply(matTransA, matL);
	matT = multiply(matInvA, matRight);

	deMat(matRight);
#endif

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
#if MASK && A_USING
	deMat(matT);
#endif
}



/*=========================================================*/
/*   MARK: Public   Functions      */
/*=========================================================*/




#if MASK
/* Encryption begins */
Mat *encrypto(
              const Mat *plain,
              const Mat *key
              
              )
{
	if (plain == NULL || key == NULL || ROUNDS < 0) return NULL;
	Mat  *matRoundKey,  *matTem;
#if A_USING
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
		sboxes(matsMasked);
		lboxes(matsMasked);
		

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

#else
/* Encryption begins */
Mat *encrypto(
        const Mat *plain,
		const Mat *key

)
{
    if(plain == NULL || key == NULL || ROUNDS < 0) return NULL;

    Mat  *roundIn, *sum, *sout, *lout;	

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


    return roundIn;
}
#endif





