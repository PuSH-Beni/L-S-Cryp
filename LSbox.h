//
//  Lbox.h
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//

#ifndef LSbox_h
#define LSbox_h


/*  - ROUNDS: The Encryption rounds
 *
 */
#define ROUNDS 16



/* -------------------------------------------------------------------------------------
 * -- MARK:   DO NOT CHANGE THE FOLLOWING SETTINGS  ------------------------------------
 * -------------------------------------------------------------------------------------
 */

#define FEISTEL 3
#include "Fundamentals.h"

/* =================================================================================
 * ============================ Constants' Defination ==============================
 * =================================================================================
 */

/* Define some matrices used in encryption, i.e. key_r, matL, rdConst etc. */
#if DIM_L == 8
#define KEY_RV  {0x68, 0x44, 0x0f, 0xc6}
#define MAT_LV  {0x89, 0x45, 0x21, 0x12, 0x78, 0x34, 0x12, 0x01 }
#define CONSTR {0}

#elif DIM_L == 16
#define KEY_RV  {	0x89, 0x45, 0x21, 0x12, 0x78, 0x34, 0x12, 0x01	}
#define MAT_LV  {	0x68, 0x44, 0x0f, 0xc6, 0xcf, 0x13, 0x9a, 0x34,\
					0x55, 0x46, 0x80, 0xa2, 0x42, 0x57, 0x09, 0x48,\
					0x77, 0xff, 0x00, 0x18, 0x27, 0x36, 0x45, 0x29,\
					0x66, 0xc7, 0x8e, 0x32, 0x47, 0x82, 0x91, 0xbb	}
#define CONSTR {	0	}

#endif //DIM_L


/* =================================================================================
 * ============================ Public Functions ===================================
 * =================================================================================
 */
#if DIM_A
Res  getMatT();
#endif
Res  encrypto(BYTE *cipher, const BYTE *plain, const BYTE *key);
Res  encrypto_fixed();
#endif /* Lbox_h */
