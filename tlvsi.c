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

#define TLVSI_CONFIG_Li			((float)3.4e-3)
#define	TLVSI_CONFIG_Lg			((float)1.8e-3)
#define TLVSI_CONFIG_Cf			((float)20e-6)

#define TLVSI_CONFIG_V_dc		((float)650.0)

#define TLVSI_CONFIG_w 			((float)314.1592653589793)

#define TLVSI_CONFIG_ts 		((float)1.0/40e3)

#define TLVSI_CONFIG_k1			((float)(TLVSI_CONFIG_ts * TLVSI_CONFIG_w))
#define TLVSI_CONFIG_k2			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Li))
#define TLVSI_CONFIG_k3			((float)((-TLVSI_CONFIG_ts / TLVSI_CONFIG_Li) * (2.0/3.0) * TLVSI_CONFIG_V_dc))
#define TLVSI_CONFIG_k4			((float)(-0.5 * TLVSI_CONFIG_ts / TLVSI_CONFIG_Cf))
#define TLVSI_CONFIG_k5			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Cf))
#define TLVSI_CONFIG_k6			((float)(-0.5 * TLVSI_CONFIG_ts / TLVSI_CONFIG_Lg))
#define TLVSI_CONFIG_k7			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Lg))

#define TLVSI_CONFIG_w_ii		((float)1.0)
#define TLVSI_CONFIG_w_ig		((float)400.0)
#define TLVSI_CONFIG_w_vc		((float)0.49)



//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
void tlvsiInitializeParams(tlvsiLCLPredictData_t *vsi, float Li, float Lg, float Cf, float V_dc, float w, float ts, float w_ii, float w_ig, float w_vc){

//	vsi->Li = Li;
//	vsi->Lg = Lg;
//	vsi->Cf = Cf;
//
//	vsi->V_dc = V_dc;
//
//	vsi->w = w;
//
//	vsi->ts = ts;
//
//	TLVSI_CONFIG_k1 = ts * w;
//	TLVSI_CONFIG_k2 = ts / Li;
//	TLVSI_CONFIG_k3 = (-ts / Li) * (2.0f / 3.0f) * V_dc;
//	TLVSI_CONFIG_k4 = -0.5f * ts / Cf;
//	TLVSI_CONFIG_k5 = ts / Cf;
//	TLVSI_CONFIG_k6 = -0.5f * ts / Lg;
//	TLVSI_CONFIG_k7 = ts / Lg;

//	vsi->w_ii = w_ii;
//	vsi->w_ig = w_ig;
//	vsi->w_vc = w_vc;

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
float tlvsiOpt(psdtypesABC_t *ii, psdtypesABC_t *ig, psdtypesABC_t *vc, psdtypesABC_t *vg){

//	static tlvsiLCLPredictData_t vsi = {.ig_d = {.d = 0.0f, .q = 0.0f; .z = 0.0f} };
	static tlvsiLCLPredictData_t vsi = {.sw = 0, .theta = 0,
			.ig_d = {.d = 0.0f, .q = 0.0f, .z = 0.0f},
			.ig_ref = {.d = -10.0f, .q = 0.0f, .z = 0.0f},
			.ii_d = {.d = 0.0f, .q = 0.0f, .z = 0.0f},
			.vc_d = {.d = 0.0f, .q = 0.0f, .z = 0.0f},
			.spll_3ph_1 = {.v_q = {0.0f, 0.0f}, .ylf = {0.0f, 0.0f}, .fo = 0.0f, .fn = 50.0f, .theta = {0.0f, 0.0f}, .delta_t = TLVSI_CONFIG_ts, .lpf_coeff = {.b0 = 166.877556f, .b1 = -166.322444f}}
	};

	float theta, phi;
    float J, Jk;
    float co, si;
    uint32_t k;

    float ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;

    /* Pre-computes sin and cos for DQ0 transforms */
    si = sinf(vsi.theta);
    co = cosf(vsi.theta);

    tptransformsABCDQ0(ii, &vsi.ii_d, si, co);
    tptransformsABCDQ0(ig, &vsi.ig_d, si, co);
    tptransformsABCDQ0(vc, &vsi.vc_d, si, co);
    tptransformsABCDQ0(vg, &vsi.vg_k, si, co);

    tppllRun(vsi.vg_k.q, &vsi.spll_3ph_1);
    vsi.theta = vsi.spll_3ph_1.theta[1];

    /* Delay compensation */
    tlvsiPredict(&vsi.ii_k, &vsi.ii_d, &vsi.ig_k, &vsi.ig_d,
                 &vsi.vc_k, &vsi.vc_d, &vsi.vg_k, vsi.theta, vsi.sw);

    /* References for filter cap. voltage and inverter current */
    vsi.vc_ref.d = ( TLVSI_CONFIG_w * TLVSI_CONFIG_Lg) * vsi.ig_ref.q + vsi.vg_k.d;
    vsi.vc_ref.q = (-TLVSI_CONFIG_w * TLVSI_CONFIG_Lg) * vsi.ig_ref.d + vsi.vg_k.q;

    vsi.ii_ref.d = vsi.ig_ref.d + ( TLVSI_CONFIG_w * TLVSI_CONFIG_Cf) * vsi.vc_ref.q;
    vsi.ii_ref.q = vsi.ig_ref.q + (-TLVSI_CONFIG_w * TLVSI_CONFIG_Cf) * vsi.vc_ref.d;

    /* Predicts for each possible switching combination */
    ii_d_constant = vsi.ii_k.d + TLVSI_CONFIG_k1 * vsi.ii_k.q + TLVSI_CONFIG_k2 * vsi.vc_k.d;
    ii_q_constant = vsi.ii_k.q - TLVSI_CONFIG_k1 * vsi.ii_k.d + TLVSI_CONFIG_k2 * vsi.vc_k.q;

    vsi.ii_k_1.d = ii_d_constant;
    vsi.ii_k_1.q = ii_q_constant;

    vc_d_constant = vsi.vc_k.d + TLVSI_CONFIG_k1 * vsi.vc_k.q + TLVSI_CONFIG_k4 * vsi.ii_k.d + TLVSI_CONFIG_k5 * vsi.ig_k.d;
    vc_q_constant = vsi.vc_k.q - TLVSI_CONFIG_k1 * vsi.vc_k.d + TLVSI_CONFIG_k4 * vsi.ii_k.q + TLVSI_CONFIG_k5 * vsi.ig_k.q;

    vsi.vc_k_1.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_d_constant);
    vsi.vc_k_1.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_q_constant);

    ig_d_constant = vsi.ig_k.d + TLVSI_CONFIG_k1 * vsi.ig_k.q + TLVSI_CONFIG_k6 * vsi.vc_k.d + TLVSI_CONFIG_k7 * vsi.vg_k.d;
    ig_q_constant = vsi.ig_k.q - TLVSI_CONFIG_k1 * vsi.ig_k.d + TLVSI_CONFIG_k6 * vsi.vc_k.q + TLVSI_CONFIG_k7 * vsi.vg_k.q;

    vsi.ig_k_1.d = ig_d_constant + TLVSI_CONFIG_k6 * (vsi.vc_k_1.d);
    vsi.ig_k_1.q = ig_q_constant + TLVSI_CONFIG_k6 * (vsi.vc_k_1.q);

//    tlvsiPredict(vsi, &vsi->ii_k_1, &vsi->ii_k, &vsi->ig_k_1, &vsi->ig_k,
//                 &vsi->vc_k_1, &vsi->vc_k, &vsi.vg_k, vsi->theta, 0);

    Jk = tlvsiCost(&vsi.ii_k_1, &vsi.ii_ref, &vsi.ig_k_1, &vsi.ig_ref,
                   &vsi.vc_k_1, &vsi.vc_ref);
    J = Jk;
    vsi.sw = 0;

    phi = 0;
    for(k = 1; k < 7; k++){

        theta = vsi.theta - phi;
        co = TLVSI_CONFIG_k3 * cosf(theta);
        si = TLVSI_CONFIG_k3 * sinf(theta);

        vsi.ii_k_1.d = ii_d_constant + co;
        vsi.ii_k_1.q = ii_q_constant - si;

        vsi.vc_k_1.d = vc_d_constant + TLVSI_CONFIG_k4 * (vsi.ii_k_1.d);
        vsi.vc_k_1.q = vc_q_constant + TLVSI_CONFIG_k4 * (vsi.ii_k_1.q);

        vsi.ig_k_1.d = ig_d_constant + TLVSI_CONFIG_k6 * (vsi.vc_k_1.d);
        vsi.ig_k_1.q = ig_q_constant + TLVSI_CONFIG_k6 * (vsi.vc_k_1.q);

//        tlvsiPredict(vsi, &vsi->ii_k_1, &vsi->ii_k, &vsi->ig_k_1, &vsi->ig_k,
//                     &vsi->vc_k_1, &vsi->vc_k, &vsi->vg_k, vsi->theta, k);

        Jk = tlvsiCost(&vsi.ii_k_1, &vsi.ii_ref,
                       &vsi.ig_k_1, &vsi.ig_ref, &vsi.vc_k_1, &vsi.vc_ref);

        if(Jk < J){
            J = Jk;
            vsi.sw = k;
        }

        phi += 1.0471975511965976f;
    }

    return J;
}
//---------------------------------------------------------------------------
void tlvsiPredict(psdtypesDQ0_t *ii_k_1, psdtypesDQ0_t *ii_k,
				  psdtypesDQ0_t *ig_k_1, psdtypesDQ0_t *ig_k,
				  psdtypesDQ0_t *vc_k_1, psdtypesDQ0_t *vc_k,
				  psdtypesDQ0_t *vg_k,
				  float theta, uint32_t sw){

	ii_k_1->d = ii_k->d + TLVSI_CONFIG_k1 * ii_k->q + TLVSI_CONFIG_k2 * vc_k->d;
	ii_k_1->q = ii_k->q - TLVSI_CONFIG_k1 * ii_k->d + TLVSI_CONFIG_k2 * vc_k->q;
	if(sw != 0){
		theta = theta - (sw-1)*1.0471975511965976f;
		ii_k_1->d +=  TLVSI_CONFIG_k3 * cosf(theta);
		ii_k_1->q += -TLVSI_CONFIG_k3 * sinf(theta);
	}

	vc_k_1->d = vc_k->d + TLVSI_CONFIG_k1 * vc_k->q + TLVSI_CONFIG_k4 * (ii_k_1->d + ii_k->d) + TLVSI_CONFIG_k5 * ig_k->d;
	vc_k_1->q = vc_k->q - TLVSI_CONFIG_k1 * vc_k->d + TLVSI_CONFIG_k4 * (ii_k_1->q + ii_k->q) + TLVSI_CONFIG_k5 * ig_k->q;

	ig_k_1->d = ig_k->d + TLVSI_CONFIG_k1 * ig_k->q + TLVSI_CONFIG_k6 * (vc_k_1->d + vc_k->d) + TLVSI_CONFIG_k7 * vg_k->d;
	ig_k_1->q = ig_k->q - TLVSI_CONFIG_k1 * ig_k->d + TLVSI_CONFIG_k6 * (vc_k_1->q + vc_k->q) + TLVSI_CONFIG_k7 * vg_k->q;
}
//---------------------------------------------------------------------------
float tlvsiCost(psdtypesDQ0_t *ii, psdtypesDQ0_t *ii_ref,
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
//---------------------------------------------------------------------------
//===========================================================================

