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

#define ROUNDS 8
#define DIM_N 64
#define DIM_S 8
#define DIM_L DIM_N/DIM_S

Mat *encryp(
	const Mat *plain,
	const Mat *key,
	const Mat *rdConst,
	const Mat *matL,
	const Mat *key_r,
	const int nRound
	);


void identMat64(
	BYTE *addr_st
	);


#endif /* Lbox_h */
