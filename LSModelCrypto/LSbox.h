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


#define ROUNDS  3
#define KEY_RV  {0x00, 0x00, 0x00, 0x00 }
#define MAT_LV  {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 } 

#define CONSTR1 {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80 } 
#define CONSTR2 {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x00 } 
#define CONSTR3 {0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00 } 


static 
Mat *rdConst[ROUNDS], *matL, *key_r;

static VECT matLV[] = MAT_LV;
static VECT keyRV[] = KEY_RV;
static VECT rdConst1V[] = CONSTR1;
static VECT rdConst2V[] = CONSTR2;
static VECT rdConst3V[] = CONSTR3;

void newPreCal();
void dePostCal();
Mat *encryp(const Mat *plain, const Mat *key);





#endif /* Lbox_h */
