//
//  Fundamentals.h
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#ifndef Fundamentals_h
#define Fundamentals_h



#include <stdlib.h>
#include <string.h>
//#include <stdio.h>
//#include <stdbool.h>
//#include <stdint.h>


#define LENG16 0
#define LENG4 0


#define DIM_L 8
#define DIM_S 8
#define ELEMS 64



#if LENG16
#define VECT WORD 
#define IDENT 0x8000
#define ONEV  0x0001
#define ZEROV 0x0000

#elif LENG4
#define VECT BYTE
#define IDENT 0x08
#define ONEV  0x01
#define ZEROV 0x00

#else
#define VECT BYTE
#define IDENT 0x80
#define ONEV  0x01
#define ZEROV 0x00
#endif

/*
typedef uint8_t             BYTE;
typedef uint16_t            WORD;
typedef uint32_t            DWORD;
typedef uint64_t            QWORD;
*/

 typedef unsigned char       BYTE;
 typedef unsigned short      WORD;
 typedef unsigned long       DWORD;
 typedef unsigned long long  QWORD;



typedef struct
{
	VECT *vect;
	int dim_row;
	int dim_col;
    char flags;	/* flags:
                *   '0': norm
                *   '1': unit matrix or indentit matrix
                *   '2': error
                */
}Mat;

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

Mat *newMat(int dim_row,  int dim_col, VECT *addr, BYTE flags);

void deMat(Mat *matrix);

/* bitand operation, matA . matB */
Mat *bitAnd(const Mat *matA, const Mat *matB);

/* add operation as same as  XOR */
Mat *add(const Mat *matA, const  Mat *matB);

/* transposition */
Mat *transpose(const Mat *matA);

/* matA x matB */
Mat *multiply(const Mat *matA, const Mat *matB);

/* catenate n mats through the r-dimension */
Mat *cat(const Mat **mats, const int n, const int r);

/* split a matrix to n parts through the r-dimension */
Mat **split (const Mat *matO, const int n, const int r);



#endif /* Fundamentals_h */




