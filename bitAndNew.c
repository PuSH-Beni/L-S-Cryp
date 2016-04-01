#include "bitAndNew.h"
extern
Mat *matAs, *matInvAs, *matTransAs;


Mat **ICA;
/*   (invA * C * A * X * A) * y
 * = y^T * A^T * ICAX^T
 * = y^T \mply ( ICAX \mply A^T )
 *
 * a row of ICAX^T : x^T \mply ICA
 * (ICAX^T composing many row vectors)
 */

/* ICA[s]:
 * I_1s*a_s1, I_1s*a_s2, ... , I_1s*a_sn
 * I_2s*a_s1, I_2s*a_s2, ... , I_2s*a_sn
 * ...
 * I_ms*a_s1, I_ms*a_s2, ... , I_ms*a_sn
 */
#define ROTATE_RIGHT(x, s, n) ((x) >> (n)) | ((x) << ((s) - (n)))

/* pre-evaluating the inv(A)CA[] array of matrices */
Mat **getMatICA(){
	int lengthOfX = DIM_L;
	// matrix A is a square matrix
	Mat **retMats = (Mat **)malloc(lengthOfX * sizeof(Mat*));
	int indexOfICA;
	Mat *invATran = transpose(matInvAs);
	for (indexOfICA = 0; indexOfICA < lengthOfX; ++indexOfICA){
		int i, k;
		BASE *ptrOfMatI, *ptrOfMatA, *ptrOfMatTem;
		int btsOfRow = bytesOfRow(matAs->dim_col);
		Mat *matTem = newMat(lengthOfX, lengthOfX, NULL, 0);

		ptrOfMatA = matAs->vect + indexOfICA * btsOfRow;
		ptrOfMatI = invATran->vect + indexOfICA * btsOfRow;
#if LENG8
		BASE check = 0x80;
#else			
		BASE check = 0x8000;
#endif
		for (i = 0; i < lengthOfX; ++i){
			ptrOfMatTem = matTem->vect + i * btsOfRow;
			ptrOfMatI += (i == LENGTH) ? 1 : 0;
			if (check & (*ptrOfMatI))
				for (k = 0; k < btsOfRow; ++k)	*(ptrOfMatTem + k) = *(ptrOfMatA + k);
			check = ROTATE_RIGHT(check, LENGTH, 1);
		}
		retMats[indexOfICA] = matTem;
	}
	return retMats;
}


void setup4BAN(){
	ICA = getMatICA();
}

/* matX or matY need to be a row vector */
Mat *bitAndNew_AA(const Mat *x, const Mat *y){
	int lengthOfX = DIM_L;
	int indexOfICA;
	Mat *tem = newMat(lengthOfX, lengthOfX, NULL, 0);
	int btsOfRow = bytesOfRow(lengthOfX);

	// get ICAX^T now
	/*for (indexOfICA = 0; indexOfICA < lengthOfX; ++indexOfICA){
		BASE *ptrOfRow = tem->vect + btsOfRow * indexOfICA;
		Mat *prod = multiply(x, ICA[indexOfICA]);
		int k;
		for (k = 0; k < btsOfRow; ++k){
		*(ptrOfRow + k) = *(prod->vect + k);
		}
		deMat(prod);
		}

		Mat *ICAX = transpose(tem);
		deMat(tem);*/

	// another methord to get ICAX:
	int row, col;
	Mat *ICAX = newMat(lengthOfX, lengthOfX, NULL, 0);
	int offset, cntsVect, k;
	BASE *ptrOfICA, *ptrOfX, *ptrOfRow;
	BASE vectTem;
	for (row = 0; row < lengthOfX; ++row){
		for (col = 0; col < lengthOfX; ++col){
			offset = col % LENGTH;
			cntsVect = col / LENGTH;
			ptrOfRow = ICAX->vect + row * btsOfRow + cntsVect;
			for (k = 0; k < btsOfRow; ++k){				
				ptrOfICA = ICA[col]->vect + row * btsOfRow + k;
				ptrOfX = x->vect + k;

				vectTem = (BASE)_mm_popcnt_u32(*ptrOfICA  & *ptrOfX) & 0x01;
				vectTem <<= (LENGTH - 1 - offset);

				*ptrOfRow ^= vectTem;
			}
		}
	}

	tem = multiply(ICAX, matTransAs);
	Mat *ret = multiply(y, tem);
	deMat(tem);
	//deMat(ICAX);

	return ret;
}
Mat *bitAndNew_EA(const Mat *x, const Mat *y){
	int lengthOfX = DIM_L;
	int i;
	Mat *ICX = newMat(lengthOfX, lengthOfX, NULL, 0);
	int btsOfRow = bytesOfRow(lengthOfX);

	// get ICX now
	for (i = 0; i < lengthOfX; ++i){
		BASE *ptrOfRow = ICX->vect + btsOfRow * i;
		BASE *ptrOfX = x->vect;
		BASE *ptrOfInvA = matInvAs->vect + btsOfRow * i;
		int k;
		for (k = 0; k < btsOfRow; ++k){
			*(ptrOfRow + k) = *(ptrOfInvA + k) & *(ptrOfX + k);
		}
	}

	Mat *tem = multiply(ICX, matTransAs);
	Mat *ret = multiply(y, tem);
	deMat(tem);
	deMat(ICX);

	return ret;
}