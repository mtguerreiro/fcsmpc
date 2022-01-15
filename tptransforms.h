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
#include "IQmathLib.h"
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

//void tptransformsABCDQ0Int(psdtypesABCint_t *abc, psdtypesDQ0int_t *dq0, fmint_t si, fmint_t co);
inline void tptransformsABCDQ0Int(psdtypesABCint_t *abc, psdtypesDQ0int_t *dq0, fmint_t si, fmint_t co){

    fmint_t alpha, beta;

    alpha   = _IQmpy( _IQ(0.66666666677f), ( abc->a - _IQmpy(_IQ(0.5f), (abc->b + abc->c)) ) );
    beta    = _IQmpy( _IQ(0.57735026913f), (abc->b - abc->c) );

    dq0->d  = _IQmpy( alpha, co) + _IQmpy(beta, si);
    dq0->q  = _IQmpy(-alpha, si) + _IQmpy(beta, co);
    dq0->z  = _IQmpy( _IQ(0.57735026913f), (abc->a + abc->b + abc->c));
}
//===========================================================================

#endif /* TPTRANSFORMS_H_ */
