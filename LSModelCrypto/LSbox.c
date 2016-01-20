//
//  LSbox.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#include "LSbox.h"

/*=========================================================*/
/*   Declarition     */
/*=========================================================*/
Mat  *rdConst[ROUNDS],  *key_r, *matL;
#if MASK
Mat *matT;
#endif
BYTE matLV[] = MAT_LV;
BYTE keyRV[] = KEY_RV;
BYTE rdConst1V[] = CONSTR1;
BYTE rdConst2V[] = CONSTR2;
BYTE rdConst3V[] = CONSTR3;

/*=========================================================*/
/*   Private   Functions      */
/*=========================================================*/

/* A 4-bit sbox */
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


#if MASK
/* Using matT, matL */
static
void lboxes(
        Mat **matsLin,
        Mat *matT,
        Mat *matL
        )
{
    Mat *matTem = matsLin[0];
    matsLin[0] = multiply(matTem,matT);
    deMat(matTem);
    int i;
    for(i = 0; i < MASKD; ++i){
        matTem = matsLin[i];
        matsLin[i] = multiply(matTem,matL);
        deMat(matTem);
    }
}

#else
/* DIM_L-bit L-box */
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



/* DIM_S-bit S-box */
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
    left = sbox4b(rows[0]);
	right = rows[1];
	
    /* Combine a bigger sbox, from 4-bit to 8-bit */
    int i;
    Mat *sum, *fout4b;
	/* Fiestel Struct */
    for(i = 0; i < 3; ++i)
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

	deMat(rows[0]);
    
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






/*=========================================================*/
/*   Public   Functions      */
/*=========================================================*/


/* Before encryption, do some pre-work to get the constant matrices */
#if MASK
void newPreCal(Mat *matAT, Mat *matAI)
#else
void newPreCal()
#endif
{
	
	key_r = newMat(DIM_S/2, DIM_L, keyRV, 0x03);
	rdConst[0] = newMat(DIM_S, DIM_L, rdConst1V, 0x03);
	rdConst[1] = newMat(DIM_S, DIM_L, rdConst2V, 0x03);
    rdConst[2] = newMat(DIM_S, DIM_L, rdConst3V, 0x03);
    matL = newMat(DIM_S, DIM_L, matLV, 0x03);
    
#if MASK    
    /* Get Matrix T  */
    Mat *matRight= multiply(matAT,matL);
    matT = multiply(matAI ,matRight);
    
    deMat(matRight);

#endif
}



/* After encryption, deconstruct those matrices */
void dePostCal()
{
	
	deMat(key_r);
	int i;
	for (i = 0; i < ROUNDS; ++i){
		deMat(rdConst[i]);
	}
#if MASK
	deMat(matT);
#endif
}


#if MASK
/* Encryption begins */
Mat *encrypto(
              const Mat *plain,
              const Mat *key,
              const Mat *matL,
              const Mat *matT
              
              )
{
	Mat *matRet = NULL;

	return matRet;
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





