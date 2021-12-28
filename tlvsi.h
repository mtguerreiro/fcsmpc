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
void tlvsiPredictInitializeParams(tlvsiLCLPredictData_t *vsi, float Li, float Lg, float Cf, float V_dc, float w, float ts);
//void tlvsiPredictInitializeDQ0Data(tlvsiLCLPredictData_t *vsi, psdtypesDQ0_t *ii_k[], psdtypesDQ0_t *ig_k[], psdtypesDQ0_t *vc_k[], psdtypesDQ0_t *vg_k);
void tlvsiPredictRun(tlvsiLCLPredictData_t *vsi, psdtypes2LLCLDQ0Data_t *lcl_k, psdtypes2LLCLDQ0Data_t *lcl_k_1, psdtypesDQ0_t vg_k, float theta, uint32_t sw);
//===========================================================================



#endif /* TLVSI_H_ */
