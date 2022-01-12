/*
 * tppll.h
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 *
 * The source code for the PLL implementation was extracted (and slightly
 * modified) from TI's DigitalPower library.
 */

#ifndef TPPLL_H_
#define TPPLL_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

#include "fixedmath.h"
//===========================================================================

typedef float float32_t;

//===========================================================================
/*------------------------------ Definitions ------------------------------*/
//===========================================================================
typedef struct{
 float32_t b1;
 float32_t b0;
} SPLL_3PH_SRF_LPF_COEFF;

typedef struct{
 float32_t v_q[2];     
 float32_t ylf[2];     
 float32_t fo;         
 float32_t fn;         
 float32_t theta[2];   
 float32_t delta_t;    
 SPLL_3PH_SRF_LPF_COEFF lpf_coeff;  
} SPLL_3PH_SRF;

typedef struct{
 fmint_t b1;
 fmint_t b0;
} SPLL_3PH_SRF_LPF_COEFF_INT;

typedef struct{
 fmint_t v_q[2];
 fmint_t ylf[2];
 fmint_t fo;
 fmint_t fn;
 fmint_t theta[2];
 fmint_t delta_t;
 SPLL_3PH_SRF_LPF_COEFF_INT lpf_coeff;
} SPLL_3PH_SRF_INT;
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
void tppllInit(float32_t grid_freq, float32_t delta_t, SPLL_3PH_SRF *spll_obj);
void tppllRun(float32_t v_q, SPLL_3PH_SRF *spll_obj);
void tppllRunInt(fmint_t v_q, SPLL_3PH_SRF_INT *spll_obj);
//===========================================================================

#endif /* TPPLL_H_ */
