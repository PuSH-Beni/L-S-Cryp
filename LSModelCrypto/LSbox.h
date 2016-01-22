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
#define ROUNDS 3
#define FEISTEL 3

/*=========================================================*/
/*   Definations About MATRICES    */
/*=========================================================*/

/* Define some matrices used in encryption, i.e. key_r, matL, rdConst etc. */
#define KEY_RV  {0x00, 0x00, 0x00, 0x00 }
#define MAT_LV  {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 }
#define CONSTR1 {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80 }
#define CONSTR2 {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x00 }
#define CONSTR3 {0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00 }




/*=========================================================*/
/*   Functions      */
/*=========================================================*/



Mat *encrypto(const Mat *plain, const Mat *key);


#endif /* Lbox_h */
