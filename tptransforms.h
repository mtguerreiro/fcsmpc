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
void tptransformsABCDQ0(psdtypesABC_t *abc, psdtypesDQ0_t *dq0, float si, float co);
void tptransformsABCDQ0Int(psdtypesABCint_t *abc, psdtypesDQ0int_t *dq0, fmint_t si, fmint_t co);

void tptransformsLineToPhase(psdtypesABCLine_t *line, psdtypesABC_t *phase);
//===========================================================================

#endif /* TPTRANSFORMS_H_ */
