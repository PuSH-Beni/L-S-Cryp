//
//  Fundamentals.h
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#ifndef Fundamentals_h
#define Fundamentals_h
/* C11 support or not */
#define C11_SUPPORT 0
/* While testing, toggle this option on */
#define TEST        1
/* Mask the text or not */
#define MASK        1

/* DIM_L can be '16' or '8' */
#define DIM_L	    16
#define DIM_S		8
#define ELEMS		(DIM_L * DIM_S)

#if MASK
/* Tree values of DIM_A
*  0:  masked without matrix A
*  4:  the length of matrix A is 4-bit
*  8:  the length of matrix A is 8-bit
*/
#define DIM_A       0
/* The dimension of the mask, i.e. x_1, x_2, ..., x_d */
#define MASKD       2

#if DIM_L == 16 && !DIM_A
#define DIVIDE_SLICES 2
#elif DIM_A
#define DIVIDE_SLICES (DIM_L / DIM_A)
#else
#define DIVIDE_SLICES 1
#endif

#endif /* MASK */

#include <stdlib.h>
#include <string.h>


/* If the compiler supports C11, 'stdbool.h' and 'stdint.h' can be used */
#if C11_SUPPORT
#include <stdbool.h>
#include <stdint.h>
#endif /* C11_SUPPORT */


/*=========================================================*/
/*   Definations About CONSTANTS    */
/*=========================================================*/


/* Use BYTE or WORD */
#define LENG16 0
#define LENG8  1

#if LENG8
#define LENGTH 8
/* 'IDENT' ---> Represents a identity vector, whose MSB is '1' */
#define IDENT 0x80
#define IDENT_MAT {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 }
#define ZERO_MAT {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }
/* The length(bits) of S-box and L-box */
/* 'ELEMS' indicates the whole length(bits) of a plain text( key, cipher etc.) */


#elif LENG16
#define LENGTH 16
#define IDENT  0x8000
#define INDENT_MAT {0x8000, 0x4000, 0x2000, 0x1000, 0x0800, 0x0400, 0x0200, 0x0100 }
#define ZERO_MAT {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }

#endif /* LENGTH */



#if C11_SUPPORT
typedef uint8_t             BYTE;
typedef uint16_t            WORD;
typedef uint32_t            DWORD;
typedef uint64_t            QWORD;
#else
typedef unsigned char       BYTE;
typedef unsigned short      WORD;
typedef unsigned long       DWORD;
typedef unsigned long long  QWORD;
#endif /* C11_SUPPORT */



/* Definations */
typedef struct
{
	BYTE *vect;
	int dim_row;
	int dim_col;
	char flags;	/* flags :
				   '0x00': normal
				   '0x01': unit matrix or identity matrix
				   '0x02': error
				   '0x03': 'vect' points to an existing array
				*/
}Mat;



/* Error types */
typedef enum
{
	RES_OK = 0,
	RES_INVALID_DIMENSION,
	RES_ALLOCATION_FAILED,
	RES_INVALID_POINTER,
	RES_INVALID_VECT,
	RES_INVALID_SPLIT,
	RES_INVALID_CAT

}Res;


#if TEST
/*=========================================================*/
/*    Private Functions      */
/*=========================================================*/
int *randOrder( );

Mat *transpose(const Mat *matO);

int bytesOfRow(int col);

BYTE shiftBit(BYTE orig, int i,  int j);

BYTE sumFromVect(BYTE Vect);

#if MASK

Mat **genRandMat( );

Mat *hatA( );

Mat *graveA( );

Mat *acuteA( );

Mat *tensorProduct(	const Mat *matX, const Mat *matY);

#endif /* MASK */

#endif /* TEST */


/*=========================================================*/
/*   Functions      */
/*=========================================================*/

/* Initialize a new Matrix instance */
Mat *newMat(int dim_row,  int dim_col, BYTE *addr, BYTE flags);

/* Deallocate a Matrix */
void deMat(Mat *matrix);

/* Deconstruct MASKD Mat instances */
void deMats(Mat **matrices, int cnt);

#if MASK
/* Encode(Mask) the plain text */
Mat **encode(const Mat *matPlain);

/* Decode(Unmask) the Secret text */
Mat *decode(Mat **matsSecret);

/* Refresh the old mask with a new mask */
Res refreshing(Mat **matsSecret);

/* Pre-calculate some matrices  */
void setup();

/* secProduct */
Mat **bitAndWithMask(const Mat **matEX, const Mat **matEY);

/* secAdd  */
Mat **addWithMask(const Mat **matEX, const  Mat **matEY);

#endif /* MASK */


/* Multiplication between two matrices */
/*  If MASK == '1'( masked ), multiply(X, Y) equals:   matX x matY  */
/*       Otherwise( unmask ), multiply(X, Y) equals:   matX x [ Transpose(matY) ]  */
Mat *multiply(const Mat *matX, const Mat *matY);

/* Simple bitand operation, as same as AND, i.e.'&' */
Mat *bitAnd(const Mat *matX, const Mat *matY);

/* Simple add operation, as same as  XOR, i.e. '^' */
Mat *add(const Mat *matX, const  Mat *matY);



/* Catenate n mats through the r-dimension */
Mat *cat(const Mat **mats,  int n,  int r);

/* Split a matrix to n parts through the r-dimension */
Mat **split (const Mat *matO,  int n,  int r);



#endif /* Fundamentals_h */





