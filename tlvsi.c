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
void tlvsiInitializeParams(tlvsiLCLPredictData_t *vsi, float Li, float Lg, float Cf, float V_dc, float w, float ts, float w_ii, float w_ig, float w_vc){

	vsi->Li = Li;
	vsi->Lg = Lg;
	vsi->Cf = Cf;

	vsi->V_dc = V_dc;

	vsi->w = w;

	vsi->ts = ts;

	vsi->k1 = ts * w;
	vsi->k2 = ts / Li;
	vsi->k3 = (-ts / Li) * (2.0f / 3.0f) * V_dc;
	vsi->k4 = -0.5f * ts / Cf;
	vsi->k5 = ts / Cf;
	vsi->k6 = -0.5f * ts / Lg;
	vsi->k7 = ts / Lg;

	vsi->w_ii = w_ii;
	vsi->w_ig = w_ig;
	vsi->w_vc = w_vc;
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
void tlvsiPredict(tlvsiLCLPredictData_t *vsi,
				  psdtypesDQ0_t *ii_k_1, psdtypesDQ0_t *ii_k,
				  psdtypesDQ0_t *ig_k_1, psdtypesDQ0_t *ig_k,
				  psdtypesDQ0_t *vc_k_1, psdtypesDQ0_t *vc_k,
				  psdtypesDQ0_t *vg_k,
				  float theta, uint32_t sw){

	ii_k_1->d = ii_k->d + vsi->k1 * ii_k->q + vsi->k2 * vc_k->d;
	ii_k_1->q = ii_k->q - vsi->k1 * ii_k->q + vsi->k2 * vc_k->q;
	if(sw != 0){
		theta = theta - (sw-1)*1.0471975511965976f;
		ii_k_1->d +=  vsi->k3 * cosf(theta);
		ii_k_1->q += -vsi->k3 * sinf(theta);
	}

	vc_k_1->d = vc_k->d + vsi->k1 * vc_k->q + vsi->k4 * (ii_k_1->d + ii_k->d) + vsi->k5 * ig_k->d;
	vc_k_1->q = vc_k->q - vsi->k1 * vc_k->d + vsi->k4 * (ii_k_1->q + ii_k->q) + vsi->k5 * ig_k->q;

	ig_k_1->d = ig_k->d + vsi->k1 * ig_k->q + vsi->k6 * (vc_k_1->d + vc_k->d) + vsi->k7 * vg_k->d;
	ig_k_1->q = ig_k->q - vsi->k1 * ig_k->d + vsi->k6 * (vc_k_1->q + vc_k->q) + vsi->k7 * vg_k->q;
}
//---------------------------------------------------------------------------
float tlvsiCost(tlvsiLCLPredictData_t *vsi,
				psdtypesDQ0_t *ii, psdtypesDQ0_t *ii_ref,
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

	J = vsi->w_ii * ((e_ii[0] * e_ii[0]) + (e_ii[1] * e_ii[1]))
	  + vsi->w_ig * ((e_ig[0] * e_ig[0]) + (e_ig[1] * e_ig[1]))
	  + vsi->w_vc * ((e_vc[0] * e_vc[0]) + (e_vc[1] * e_vc[1]));

	return J;
}
//---------------------------------------------------------------------------
//===========================================================================

