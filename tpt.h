/*
 * tptransforms.h
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 *
 * Three-phase transforms.
 *
 * The source code for the transform was extracted (and slightly modified)
 * from TI's DigitalPower library.
 */

#ifndef TPT_H_
#define TPT_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

#include "tpdt.h"
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
void tptAbcDq0(tpdtAbc_t *abc, tpdtDq0_t *dq0, tpdtfloat_t si, tpdtfloat_t co);
void tptAbcDq0FixedPoint(tpdtAbcFixedPoint_t *abc, tpdtDq0FixedPoint_t *dq0, tpdtint_t si, tpdtint_t co);

void tptLineToPhase(tpdtAbcLine_t *line, tpdtAbc_t *phase);
void tptLineToPhaseFixedPoint(tpdtAbcLineFixedPoint_t *line, tpdtAbcFixedPoint_t *phase);
//===========================================================================

#endif /* TPT_H_ */
