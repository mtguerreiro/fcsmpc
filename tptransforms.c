/*
 * tptransforms.c
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 */

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include "tptransforms.h"

#include <math.h>
//===========================================================================


#define TPTRANSFORMS_C1_INT		((int32_t)(0.66666666677f * 65536.0f))
#define TPTRANSFORMS_C2_INT		((int32_t)(0.5f * 65536.0f))
#define TPTRANSFORMS_C3_INT		((int32_t)(0.57735026913f * 65536.0f))

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
void tptransformsABCDQ0(psdtypesABC_t *abc, psdtypesDQ0_t *dq0, float si, float co){

	float alpha, beta;

    alpha	= (0.66666666677f) * (abc->a - 0.5f * (abc->b + abc->c));
    beta	= (0.57735026913f) * (abc->b - abc->c);

    dq0->d	=  alpha * co + beta * si;
    dq0->q	= -alpha * si + beta * co;
    dq0->z	= (0.57735026913f) * (abc->a + abc->b + abc->c);
}
//---------------------------------------------------------------------------
void tptransformsABCDQ0Fixed(psdtypesABCint_t *abc, psdtypesDQ0int_t *dq0, int32_t si, int32_t co){

	int32_t alpha, beta;

    alpha	= (TPTRANSFORMS_C1_INT * (abc->a - ((TPTRANSFORMS_C2_INT * (abc->b + abc->c)) >> 16))) >> 16;
    beta	= ((TPTRANSFORMS_C3_INT) * (abc->b - abc->c)) >> 16;

    dq0->d	=  ((alpha * co) >> 16) + ((beta * si) >> 16);
    dq0->q	= ((-alpha * si) >> 16) + ((beta * co) >> 16);
    dq0->z	= ((TPTRANSFORMS_C3_INT) * (abc->a + abc->b + abc->c)) >> 16;
}
//===========================================================================

