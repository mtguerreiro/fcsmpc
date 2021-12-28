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
//===========================================================================

//===========================================================================
/*------------------------------ Definitions ------------------------------*/
//===========================================================================
typedef struct{
	float ts;

	float Li;
	float Lg;
	float Cf;

	float w;

	float V_dc;

	float k1, k2, k3, k4, k5, k6, k7;

	float w_ii, w_ig, w_vc;

//	psdtypesDQ0_t *ii_k[2];
//
//	psdtypesDQ0_t *ig_k[2];
//
//	psdtypesDQ0_t *vc_k[2];
//
//	psdtypesDQ0_t *vg_k;
} tlvsiLCLPredictData_t;
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
void tlvsiInitializeParams(tlvsiLCLPredictData_t *vsi, float Li, float Lg, float Cf, float V_dc, float w, float ts, float w_ii, float w_ig, float w_vc);
//void tlvsiPredictInitializeDQ0Data(tlvsiLCLPredictData_t *vsi, psdtypesDQ0_t *ii_k[], psdtypesDQ0_t *ig_k[], psdtypesDQ0_t *vc_k[], psdtypesDQ0_t *vg_k);
void tlvsiPredict(tlvsiLCLPredictData_t *vsi,
				  psdtypesDQ0_t *ii_k_1, psdtypesDQ0_t *ii_k,
				  psdtypesDQ0_t *ig_k_1, psdtypesDQ0_t *ig_k,
				  psdtypesDQ0_t *vc_k_1, psdtypesDQ0_t *vc_k,
				  psdtypesDQ0_t *vg_k,
				  float theta, uint32_t sw);

float tlvsiCost(tlvsiLCLPredictData_t *vsi,
				psdtypesDQ0_t *ii, psdtypesDQ0_t *ii_ref,
				psdtypesDQ0_t *ig, psdtypesDQ0_t *ig_ref,
				psdtypesDQ0_t *vc, psdtypesDQ0_t *vc_ref);
//===========================================================================



#endif /* TLVSI_H_ */
