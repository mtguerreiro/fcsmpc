/*
 * tpt.c
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 */

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include "tpt.h"

#include <math.h>

#include "fixedmath.h"
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
void tptAbcDq0(tpdtAbc_t *abc, tpdtDq0_t *dq0, tpdtfloat_t si, tpdtfloat_t co){

	float alpha, beta;

    alpha	= (0.66666666677f) * (abc->a - 0.5f * (abc->b + abc->c));
    beta	= (0.57735026913f) * (abc->b - abc->c);

    dq0->d	=  alpha * co + beta * si;
    dq0->q	= -alpha * si + beta * co;
    dq0->z	= (0.57735026913f) * (abc->a + abc->b + abc->c);
}
//---------------------------------------------------------------------------
void tptAbcDq0FixedPoint(tpdtAbcFixedPoint_t *abc, tpdtDq0FixedPoint_t *dq0, tpdtint_t si, tpdtint_t co){

	fmint_t alpha, beta;

    alpha	= fixedmul( fixedmathftoi(0.66666666677f), ( abc->a - fixedmul(fixedmathftoi(0.5f), (abc->b + abc->c)) ) );
    beta	= fixedmul( fixedmathftoi(0.57735026913f), (abc->b - abc->c) );

    dq0->d	= fixedmul( alpha, co) + fixedmul(beta, si);
    dq0->q	= fixedmul(-alpha, si) + fixedmul(beta, co);
    dq0->z	= fixedmul( fixedmathftoi(0.57735026913f), (abc->a + abc->b + abc->c));
}
//---------------------------------------------------------------------------
void tptLineToPhase(tpdtAbcLine_t *line, tpdtAbc_t *phase){

    phase->a = (line->ac + line->ab) / 3.0f;
    phase->b = (line->ac - 2.0f * line->ab) / 3.0f;
    phase->c = (-2.0f * line->ac + line->ab) / 3.0f;
}
//---------------------------------------------------------------------------
void tptLineToPhaseFixedPoint(tpdtAbcLineFixedPoint_t *line, tpdtAbcFixedPoint_t *phase){

    phase->a = fixedmul( (line->ac + line->ab), fixedmathftoi(1.0f / 3.0f) );
    phase->b = fixedmul( (line->ac - fixedmul(fixedmathftoi(2.0f), line->ab)), fixedmathftoi(1.0f / 3.0f) );
    phase->c = fixedmul( (-fixedmul(fixedmathftoi(2.0f), line->ac) + line->ab), fixedmathftoi(1.0f / 3.0f) );
}
//---------------------------------------------------------------------------
//===========================================================================

