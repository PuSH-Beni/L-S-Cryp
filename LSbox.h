//
//  Lbox.h
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#ifndef LSbox_h
#define LSbox_h

#include "Fundamentals.h"

/*=========================================================*/
/*   Definations     */
/*=========================================================*/
#define ROUNDS 16
#define FEISTEL 3

/* About MATRICES    */
/* Define some matrices used in encryption, i.e. key_r, matL, rdConst etc. */
#if DIM_L == 8
#define KEY_RV  {0x68, 0x44, 0x0f, 0xc6}
#define MAT_LV  {0x89, 0x45, 0x21, 0x12, 0x78, 0x34, 0x12, 0x01 }

#define CONSTR1 {0}
#define CONSTR2 {0}
#define CONSTR3 {0}
#define CONSTR4 {0}
#define CONSTR5 {0}
#define CONSTR6 {0}
#define CONSTR7 {0}
#define CONSTR8 {0}
#define CONSTR9 {0}
#define CONSTR10 {0}
#define CONSTR11 {0}
#define CONSTR12 {0}
#define CONSTR13 {0}
#define CONSTR14 {0}
#define CONSTR15 {0}
#define CONSTR16 {0}

#elif DIM_L == 16
#define KEY_RV  {0x89, 0x45, 0x21, 0x12, 0x78, 0x34, 0x12, 0x01}
#define MAT_LV  {0x68, 0x44, 0x0f, 0xc6,  0xcf, 0x13, 0x9a, 0x34, 0x55, 0x46, 0x80, 0xa2, 0x42, 0x57, 0x09, 0x48,0x77, 0xff, 0x00, 0x18, 0x27, 0x36, 0x45, 0x29, 0x66, 0xc7, 0x8e, 0x32, 0x47, 0x82, 0x91, 0xbb}

#define CONSTR1 {0}
#define CONSTR2 {0}
#define CONSTR3 {0}
#define CONSTR4 {0}
#define CONSTR5 {0}
#define CONSTR6 {0}
#define CONSTR7 {0}
#define CONSTR8 {0}
#define CONSTR9  {0}
#define CONSTR10 {0}
#define CONSTR11 {0}
#define CONSTR12 {0}
#define CONSTR13 {0}
#define CONSTR14 {0}
#define CONSTR15 {0}
#define CONSTR16 {0}
#endif


#define KEY_SIZE (DIM_S / 2 * (DIM_L / 8))
#define L_SIZE   (DIM_L * (DIM_L / 8))
#define CONST_SIZE  (DIM_S * (DIM_L / 8))
/*=========================================================*/
/*   Functions      */
/*=========================================================*/
void enPreCal();
void dePostCal();
Mat *encrypto(const Mat *plain, const Mat *key);

#endif /* Lbox_h */
