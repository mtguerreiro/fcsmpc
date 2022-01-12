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

    float theta;

    uint32_t sw;

    SPLL_3PH_SRF spll_3ph_1;

} tlvsiLCLPredictDataInt_t;
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

uint32_t tlvsiOptFixed(psdtypesDQ0int_t *ii, psdtypesDQ0int_t *ig, psdtypesDQ0int_t *vc, psdtypesDQ0int_t *vg, float theta, fmint_t *Jopt);

float tlvsiCost(psdtypesDQ0_t *ii, psdtypesDQ0_t *ii_ref,
				psdtypesDQ0_t *ig, psdtypesDQ0_t *ig_ref,
				psdtypesDQ0_t *vc, psdtypesDQ0_t *vc_ref);

fmint_t tlvsiCostFixed(psdtypesDQ0int_t *ii, psdtypesDQ0int_t *ii_ref,
					 psdtypesDQ0int_t *ig, psdtypesDQ0int_t *ig_ref,
					 psdtypesDQ0int_t *vc, psdtypesDQ0int_t *vc_ref);
//===========================================================================



#endif /* TLVSI_H_ */
