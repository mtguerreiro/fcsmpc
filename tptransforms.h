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
#include <psdtypes.h>
#include <stdint.h>

//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
void tptransformsABCDQ0(psdtypesABC_t *abc, psdtypesDQ0_t *dq0, float theta);
//===========================================================================

#endif /* TPTRANSFORMS_H_ */
