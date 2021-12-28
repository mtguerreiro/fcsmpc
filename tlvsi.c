/*
 * tlvsi.c
 *
 *  Created on: 28 de dez de 2021
 *      Author: marco
 */

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include "tlvsi.h"

#include <math.h>
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
void tlvsiPredictInitializeParams(tlvsiLCLPredictData_t *vsi, float Li, float Lg, float Cf, float V_dc, float w, float ts){

	vsi->Li = Li;
	vsi->Lg = Lg;
	vsi->Cf = Cf;

	vsi->V_dc = V_dc;

	vsi->w = w;

	vsi->ts = ts;

	vsi->k1 = ts * w;
	vsi->k2 = ts / Li;
	vsi->k3 = (-ts / Li) * (2.0f / 3.0f) * V_dc;
	vsi->k4 = 0.5f * ts / Cf;
	vsi->k5 = ts / Cf;
	vsi->k6 = -0.5f * ts / Lg;
	vsi->k7 = ts / Lg;
}
//---------------------------------------------------------------------------
//void tlvsiPredictInitializeDQ0Data(tlvsiLCLPredictData_t *vsi, psdtypesDQ0_t *ii_k[], psdtypesDQ0_t *ig_k[], psdtypesDQ0_t *vc_k[], psdtypesDQ0_t *vg_k){
//
//	vsi->ii_k = ii_k;
//	vsi->ig_k = ig_k;
//	vsi->vc_k = vc_k;
//	vsi->vg_k = vg_k;
//}
//---------------------------------------------------------------------------
void tlvsiPredictRun(tlvsiLCLPredictData_t *vsi, psdtypes2LLCLDQ0Data_t *lcl_k, psdtypes2LLCLDQ0Data_t *lcl_k_1, psdtypesDQ0_t vg_k, float theta, uint32_t sw){


	lcl_k_1->ii.d = lcl_k->ii.d + vsi->k1 * lcl_k->ii.q + vsi->k2 * lcl_k->vc.d;
	lcl_k_1->ii.q = lcl_k->ii.q - vsi->k1 * lcl_k->ii.q + vsi->k2 * lcl_k->vc.q;
	if(sw != 0){
		theta = theta - (sw-1)*1.0471975511965976f;
		lcl_k_1->ii.d +=  vsi->k3 * cosf(theta);
		lcl_k_1->ii.q += -vsi->k3 * sinf(theta);
	}

	lcl_k_1->vc.d = lcl_k->vc.d + vsi->k1 * lcl_k->vc.q + vsi->k4 * (lcl_k_1->ii.d + lcl_k->ii.d) + vsi->k5 * lcl_k->ig.d;
	lcl_k_1->vc.q = lcl_k->vc.q - vsi->k1 * lcl_k->vc.d + vsi->k4 * (lcl_k_1->ii.q + lcl_k->ii.q) + vsi->k5 * lcl_k->ig.q;

	lcl_k_1->ig.d = lcl_k->ig.d + vsi->k1 * lcl_k->ig.q + vsi->k6 * (lcl_k_1->vc.d + lcl_k->vc.d) + vsi->k7 * vg_k.d;
	lcl_k_1->ig.q = lcl_k->ig.q - vsi->k1 * lcl_k->ig.d + vsi->k4 * (lcl_k_1->vc.q + lcl_k->vc.q) + vsi->k7 * vg_k.q;
}
//---------------------------------------------------------------------------
//===========================================================================

