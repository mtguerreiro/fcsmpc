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

#include "tptransforms.h"
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

    vsi->ig_d.d = 0.0f;
    vsi->ig_d.q = 0.0f;
    vsi->ig_d.z = 0.0f;

    vsi->ig_ref.d = -10.0f;
    vsi->ig_ref.q = 0.0f;
    vsi->ig_ref.z = 0.0f;

    vsi->ii_d.d = 0.0f;
    vsi->ii_d.q = 0.0f;
    vsi->ii_d.z = 0.0f;

    vsi->vc_d.d = 0.0f;
    vsi->vc_d.q = 0.0f;
    vsi->vc_d.z = 0.0f;

    vsi->theta = 0.0f;
    vsi->sw = 0;

    tppllInit(50.0f, ts, &vsi->spll_3ph_1);
    vsi->spll_3ph_1.lpf_coeff.b0 = 166.877556f;
    vsi->spll_3ph_1.lpf_coeff.b1 = -166.322444f;
}
//---------------------------------------------------------------------------
float tlvsiOpt(tlvsiLCLPredictData_t *vsi,
               psdtypesABC_t *ii, psdtypesABC_t *ig,
               psdtypesABC_t *vc, psdtypesABC_t *vg){

    float J, Jk;
    uint32_t k;

    tptransformsABCDQ0(ii, &vsi->ii_d, vsi->theta);
    tptransformsABCDQ0(ig, &vsi->ig_d, vsi->theta);
    tptransformsABCDQ0(vc, &vsi->vc_d, vsi->theta);
    tptransformsABCDQ0(vg, &vsi->vg_k, vsi->theta);

    tppllRun(vsi->vg_k.q, &vsi->spll_3ph_1);
    vsi->theta = vsi->spll_3ph_1.theta[1];

    /* Delay compensation */
    tlvsiPredict(vsi, &vsi->ii_k, &vsi->ii_d, &vsi->ig_k, &vsi->ig_d,
                 &vsi->vc_k, &vsi->vc_d, &vsi->vg_k, vsi->theta, vsi->sw);

    /* References for filter cap. voltage and inverter current */
    vsi->vc_ref.d = ( vsi->w * vsi->Lg) * vsi->ig_ref.q + vsi->vg_k.d;
    vsi->vc_ref.q = (-vsi->w * vsi->Lg) * vsi->ig_ref.d + vsi->vg_k.q;

    vsi->ii_ref.d = vsi->ig_ref.d + ( vsi->w * vsi->Cf) * vsi->vc_ref.q;
    vsi->ii_ref.q = vsi->ig_ref.q + (-vsi->w * vsi->Cf) * vsi->vc_ref.d;

    /* Predicts for each possible switching combination */
    tlvsiPredict(vsi, &vsi->ii_k_1, &vsi->ii_k, &vsi->ig_k_1, &vsi->ig_k,
                 &vsi->vc_k_1, &vsi->vc_k, &vsi->vg_k, vsi->theta, 0);
    Jk = tlvsiCost(vsi, &vsi->ii_k_1, &vsi->ii_ref, &vsi->ig_k_1, &vsi->ig_ref,
                   &vsi->vc_k_1, &vsi->vc_ref);
    J = Jk;
    vsi->sw = 0;

    for(k = 1; k < 7; k++){

        tlvsiPredict(vsi, &vsi->ii_k_1, &vsi->ii_k, &vsi->ig_k_1, &vsi->ig_k,
                     &vsi->vc_k_1, &vsi->vc_k, &vsi->vg_k, vsi->theta, k);

        Jk = tlvsiCost(vsi, &vsi->ii_k_1, &vsi->ii_ref,
                       &vsi->ig_k_1, &vsi->ig_ref, &vsi->vc_k_1, &vsi->vc_ref);

        if(Jk < J){
            J = Jk;
            vsi->sw = k;
        }
    }

    return J;
}
//---------------------------------------------------------------------------
void tlvsiPredict(tlvsiLCLPredictData_t *vsi,
				  psdtypesDQ0_t *ii_k_1, psdtypesDQ0_t *ii_k,
				  psdtypesDQ0_t *ig_k_1, psdtypesDQ0_t *ig_k,
				  psdtypesDQ0_t *vc_k_1, psdtypesDQ0_t *vc_k,
				  psdtypesDQ0_t *vg_k,
				  float theta, uint32_t sw){

	ii_k_1->d = ii_k->d + vsi->k1 * ii_k->q + vsi->k2 * vc_k->d;
	ii_k_1->q = ii_k->q - vsi->k1 * ii_k->d + vsi->k2 * vc_k->q;
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

