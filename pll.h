/*
 * pll.h
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 *
 * The source code for the PLL implementation was extracted (and slightly
 * modified) from TI's DigitalPower library.
 */

#ifndef PLL_H_
#define PLL_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

#include "tpdt.h"
//===========================================================================

typedef tpdtfloat_t pllfloat_t;
typedef tpdtint_t pllint_t;

//===========================================================================
/*------------------------------ Definitions ------------------------------*/
//===========================================================================
typedef struct{
	pllfloat_t b1;
	pllfloat_t b0;
}pllLpfCoeff;

typedef struct{
	pllfloat_t v_q[2];
	pllfloat_t ylf[2];
	pllfloat_t fo;
	pllfloat_t fn;
	pllfloat_t theta[2];
	pllfloat_t delta_t;
	pllLpfCoeff lpf_coeff;
}pll_t;

typedef struct{
	pllint_t b1;
	pllint_t b0;
}pllLpfCoeffInt;

typedef struct{
	pllint_t v_q[2];
	pllint_t ylf[2];
	pllint_t fo;
	pllint_t fn;
	pllint_t theta[2];
	pllint_t delta_t;
	pllLpfCoeffInt lpf_coeff;
}pllFixedPoint_t;
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
void pllInit(pll_t *pll, pllfloat_t fn, pllfloat_t dt);
void pllRun(pll_t *pll, pllfloat_t v_q);

void pllInitFixedPoint(pllFixedPoint_t *pll, pllint_t fn, pllint_t dt);
void pllRunFixedPoint(pllFixedPoint_t *pll, pllint_t v_q);
//===========================================================================

#endif /* PLL_H_ */
