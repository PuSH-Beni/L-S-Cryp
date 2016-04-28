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
 *
 */

#define DIM_A       8


/* The dimension of the mask, i.e. x_1, x_2, ..., x_d */
#define MASKD       2
#define MASKD_SQURE MASKD * MASKD
#define DIVIDE      0

#if DIM_A
#define DIVIDE_PARTS (DIM_L / DIM_A)
#else
#define DIVIDE_PARTS 1
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
#define LENG8  1

#if LENG8
#define LENGTH 8
/* 'IDENT' ---> Represents a identity vector, whose MSB is '1' */
#define UNIT_BYTE 0x80
#define UNIT_MAT {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 }
#define ZERO_MAT  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }
/* The length(bits) of S-box and L-box */
/* 'ELEMS' indicates the whole length(bits) of a plain text( key, cipher etc.) */

#elif LENG16

#define LENGTH 16
#define IDENT  0x8000
#define INDENT_MAT {0x8000, 0x4000, 0x2000, 0x1000, 0x0800, 0x0400, 0x0200, 0x0100 }
#define ZERO_MAT {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }
typedef unsigned short BASE;
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



/* Error types */
typedef enum
{
	RES_OK = 0,
	RES_INVALID_DIMENSION,
	RES_INVALID_POINTER,
	RES_INVALID_VECT,
	RES_NON_MATA

}Res;





/*=========================================================*/
/*   Functions      */
/*=========================================================*/

int bytesOfRow(int col);
#if MASK
/* Encode(Mask) the plain text */
Res encode(BYTE *matMasked, const BYTE *matPlain);

/* Decode(Unmask) the Secret text */
Res decode(BYTE *matUnmask, const BYTE *matsSecret);

/* Pre-calculate some matrices  */
Res setupEnc();

/* secProduct */
Res bitAndWithMask(BYTE *bitAndRes, const BYTE *matEX, const BYTE *matEY, const int *dims);

/* secAdd  */
Res addWithMask(BYTE *addRes, const BYTE *matEX, const  BYTE *matEY, const int *dims);

#endif /* MASK */


/* Multiply(X, Y) equals:   matX x [ Transpose(matY) ]  */
Res multiply(BYTE *multiRes, const BYTE *matX, const BYTE *matY, const int *dims);

/* Simple bitand operation, as same as AND, i.e.'&' */
Res bitAnd(BYTE *bitAndRes, const BYTE *matX, const BYTE *matY, const int *dims);

/* Simple add operation, as same as  XOR, i.e. '^' */
Res add(BYTE *addRes, const BYTE *matX, const  BYTE *matY, const int *dims);

Res transpose(BYTE *transRes, const BYTE *matOrig, const int *dims);


#endif /* Fundamentals_h */





