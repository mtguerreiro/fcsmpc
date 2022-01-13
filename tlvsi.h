/*
 * tlvsi.h
 *
 *  Created on: 28 de dez de 2021
 *      Author: marco
 */

#ifndef TLVSI_H_
#define TLVSI_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

#include "psdtypes.h"
#include "tppll.h"

#include "fixedmath.h"
//===========================================================================

//===========================================================================
/*------------------------------ Definitions ------------------------------*/
//===========================================================================
typedef struct{

    psdtypesABC_t ig_abc;
    psdtypesDQ0_t ig_d;
    psdtypesDQ0_t ig_k;
    psdtypesDQ0_t ig_k_1;
    psdtypesDQ0_t ig_ref;

    psdtypesABC_t ii_abc;
    psdtypesDQ0_t ii_d;
    psdtypesDQ0_t ii_k;
    psdtypesDQ0_t ii_k_1;
    psdtypesDQ0_t ii_ref;

    psdtypesABC_t vc_abc;
    psdtypesDQ0_t vc_d;
    psdtypesDQ0_t vc_k;
    psdtypesDQ0_t vc_k_1;
    psdtypesDQ0_t vc_ref;

    psdtypesABC_t vg_abc;
    psdtypesDQ0_t vg_k;

    float theta;

    uint32_t sw;

    SPLL_3PH_SRF spll_3ph_1;

} tlvsiLCLPredictData_t;

typedef struct{

    psdtypesABC_t ig_abc;
    psdtypesDQ0int_t ig_d;
    psdtypesDQ0int_t ig_k;
    psdtypesDQ0int_t ig_k_1;
    psdtypesDQ0int_t ig_ref;

    psdtypesABC_t ii_abc;
    psdtypesDQ0int_t ii_d;
    psdtypesDQ0int_t ii_k;
    psdtypesDQ0int_t ii_k_1;
    psdtypesDQ0int_t ii_ref;

    psdtypesABC_t vc_abc;
    psdtypesDQ0int_t vc_d;
    psdtypesDQ0int_t vc_k;
    psdtypesDQ0int_t vc_k_1;
    psdtypesDQ0int_t vc_ref;

    psdtypesABC_t vg_abc;
    psdtypesDQ0int_t vg_k;

    fmint_t theta;

    uint32_t sw;

    SPLL_3PH_SRF_INT spll_3ph_1;

} tlvsiLCLPredictDataInt_t;

#define TLVSI_CONFIG_w_ii		((float)(1.0f / 100.0f))
#define TLVSI_CONFIG_w_ig		((float)(400.0f / 100.0f))
#define TLVSI_CONFIG_w_vc		((float)(0.49f / 100.0f))

#define TLVSI_CONFIG_w_ii_int	fixedmathftoi(TLVSI_CONFIG_w_ii)
#define TLVSI_CONFIG_w_ig_int	fixedmathftoi(TLVSI_CONFIG_w_ig)
#define TLVSI_CONFIG_w_vc_int	fixedmathftoi(TLVSI_CONFIG_w_vc)
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
uint32_t tlvsiOpt(psdtypesABC_t *ii, psdtypesABC_t *ig, psdtypesABC_t *vc, psdtypesABC_t *vg, float *Jopt);

void tlvsiPredict(psdtypesDQ0_t *ii_k_1, psdtypesDQ0_t *ii_k,
				  psdtypesDQ0_t *ig_k_1, psdtypesDQ0_t *ig_k,
				  psdtypesDQ0_t *vc_k_1, psdtypesDQ0_t *vc_k,
				  psdtypesDQ0_t *vg_k,
				  float theta, uint32_t sw);

void tlvsiPredictFixed(psdtypesDQ0int_t *ii_k_1, psdtypesDQ0int_t *ii_k,
				  	   psdtypesDQ0int_t *ig_k_1, psdtypesDQ0int_t *ig_k,
					   psdtypesDQ0int_t *vc_k_1, psdtypesDQ0int_t *vc_k,
					   psdtypesDQ0int_t *vg_k,
					   fmint_t theta, uint32_t sw);

uint32_t tlvsiOptFixed(psdtypesABCint_t *ii, psdtypesABCint_t *ig, psdtypesABCint_t *vc, psdtypesABCint_t *vg, fmint_t *Jopt);

//float tlvsiCost(psdtypesDQ0_t *ii, psdtypesDQ0_t *ii_ref,
//				psdtypesDQ0_t *ig, psdtypesDQ0_t *ig_ref,
//				psdtypesDQ0_t *vc, psdtypesDQ0_t *vc_ref);
inline float tlvsiCost(psdtypesDQ0_t *ii, psdtypesDQ0_t *ii_ref,
				psdtypesDQ0_t *ig, psdtypesDQ0_t *ig_ref,
				psdtypesDQ0_t *vc, psdtypesDQ0_t *vc_ref){

	float e_ii[2], e_ig[2], e_vc[2];
	float J;

	e_ii[0] = ii->d - ii_ref->d;
	e_ii[1] = ii->q - ii_ref->q;

	e_ig[0] = ig->d - ig_ref->d;
	e_ig[1] = ig->q - ig_ref->q;

	e_vc[0] = vc->d - vc_ref->d;
	e_vc[1] = vc->q - vc_ref->q;

	J = TLVSI_CONFIG_w_ii * ((e_ii[0] * e_ii[0]) + (e_ii[1] * e_ii[1]))
	  + TLVSI_CONFIG_w_ig * ((e_ig[0] * e_ig[0]) + (e_ig[1] * e_ig[1]))
	  + TLVSI_CONFIG_w_vc * ((e_vc[0] * e_vc[0]) + (e_vc[1] * e_vc[1]));

	return J;
}

fmint_t tlvsiCostFixed(psdtypesDQ0int_t *ii, psdtypesDQ0int_t *ii_ref,
					 psdtypesDQ0int_t *ig, psdtypesDQ0int_t *ig_ref,
					 psdtypesDQ0int_t *vc, psdtypesDQ0int_t *vc_ref);

//===========================================================================



#endif /* TLVSI_H_ */
