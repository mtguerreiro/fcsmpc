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
void tptransformsABCDQ0Int(psdtypesABCint_t *abc, psdtypesDQ0int_t *dq0, fmint_t si, fmint_t co){

	fmint_t alpha, beta;

    alpha	= fixedmul( fixedmathftoi(0.66666666677f), ( abc->a - fixedmul(fixedmathftoi(0.5f), (abc->b + abc->c)) ) );
    beta	= fixedmul( fixedmathftoi(0.57735026913f), (abc->b - abc->c) );

    dq0->d	= fixedmul( alpha, co) + fixedmul(beta, si);
    dq0->q	= fixedmul(-alpha, si) + fixedmul(beta, co);
    dq0->z	= fixedmul( fixedmathftoi(0.57735026913f), (abc->a + abc->b + abc->c));
}
//===========================================================================

