/*
 * tptransforms.h
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 *
 * The source code for the transform was extracted (and slightly modified)
 * from TI's DigitalPower library.
 */

#ifndef TPTRANSFORMS_H_
#define TPTRANSFORMS_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

#include "psdtypes.h"

#include "fixedmath.h"
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//void tptransformsABCDQ0(psdtypesABC_t *abc, psdtypesDQ0_t *dq0, float si, float co);
inline void tptransformsABCDQ0(psdtypesABC_t *abc, psdtypesDQ0_t *dq0, float si, float co){

    float alpha, beta;

    alpha   = (0.66666666677f) * (abc->a - 0.5f * (abc->b + abc->c));
    beta    = (0.57735026913f) * (abc->b - abc->c);

    dq0->d  =  alpha * co + beta * si;
    dq0->q  = -alpha * si + beta * co;
    dq0->z  = (0.57735026913f) * (abc->a + abc->b + abc->c);
}

void tptransformsABCDQ0Int(psdtypesABCint_t *abc, psdtypesDQ0int_t *dq0, fmint_t si, fmint_t co);
//===========================================================================

#endif /* TPTRANSFORMS_H_ */
