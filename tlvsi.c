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

#define TLVSI_CONFIG_V_dc		((float)(650.0f / 200.0f))

#define TLVSI_CONFIG_w 			((float)314.1592653589793)

#define TLVSI_CONFIG_ts 		((float)1.0f/40e3)

#define TLVSI_CONFIG_k1			((float)(TLVSI_CONFIG_ts * TLVSI_CONFIG_w))
#define TLVSI_CONFIG_k2			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Li))
#define TLVSI_CONFIG_k3			((float)((-TLVSI_CONFIG_ts / TLVSI_CONFIG_Li) * (2.0f/3.0f) * TLVSI_CONFIG_V_dc))
#define TLVSI_CONFIG_k4			((float)(-0.5f * TLVSI_CONFIG_ts / TLVSI_CONFIG_Cf))
#define TLVSI_CONFIG_k5			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Cf))
#define TLVSI_CONFIG_k6			((float)(-0.5f * TLVSI_CONFIG_ts / TLVSI_CONFIG_Lg))
#define TLVSI_CONFIG_k7			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Lg))

#define TLVSI_CONFIG_k8			(TLVSI_CONFIG_w * TLVSI_CONFIG_Lg)
#define TLVSI_CONFIG_k9			(TLVSI_CONFIG_w * TLVSI_CONFIG_Cf)

#define TLVSI_CONFIG_w_ii		((float)(1.0f / 100.0f))
#define TLVSI_CONFIG_w_ig		((float)(400.0f / 100.0f))
#define TLVSI_CONFIG_w_vc		((float)(0.49f / 100.0f))

#define TLVSI_CONFIG_k1_int		fixedmathftoi(TLVSI_CONFIG_k1)
#define TLVSI_CONFIG_k2_int		fixedmathftoi(TLVSI_CONFIG_k2)
#define TLVSI_CONFIG_k3_int		fixedmathftoi(TLVSI_CONFIG_k3)
#define TLVSI_CONFIG_k4_int		fixedmathftoi(TLVSI_CONFIG_k4)
#define TLVSI_CONFIG_k5_int		fixedmathftoi(TLVSI_CONFIG_k5)
#define TLVSI_CONFIG_k6_int		fixedmathftoi(TLVSI_CONFIG_k6)
#define TLVSI_CONFIG_k7_int		fixedmathftoi(TLVSI_CONFIG_k7)
#define TLVSI_CONFIG_k8_int		fixedmathftoi(TLVSI_CONFIG_k8)
#define TLVSI_CONFIG_k9_int		fixedmathftoi(TLVSI_CONFIG_k9)

#define TLVSI_CONFIG_w_ii_int	fixedmathftoi(TLVSI_CONFIG_w_ii)
#define TLVSI_CONFIG_w_ig_int	fixedmathftoi(TLVSI_CONFIG_w_ig)
#define TLVSI_CONFIG_w_vc_int	fixedmathftoi(TLVSI_CONFIG_w_vc)

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
uint32_t tlvsiOpt(psdtypesABC_t *ii, psdtypesABC_t *ig, psdtypesABC_t *vc, psdtypesABC_t *vg, float *Jopt){

	static tlvsiLCLPredictData_t vsi = {.sw = 0, .theta = 0.0f,
			.ig_d = {.d = 0.0f, .q = 0.0f, .z = 0.0f},
			.ig_ref = {.d = -10.0f / 200.0f, .q = 0.0f, .z = 0.0f},
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

    if( Jopt != 0 ) *Jopt = J;

    return vsi.sw;
}
//---------------------------------------------------------------------------
uint32_t tlvsiOptFixed(psdtypesDQ0int_t *ii, psdtypesDQ0int_t *ig, psdtypesDQ0int_t *vc, psdtypesDQ0int_t *vg, float theta, fmint_t *Jopt){

	static tlvsiLCLPredictDataInt_t vsi = {.sw = 0, .theta = 0.0f,
			.ig_d = {.d = 0, .q = 0, .z = 0},
			.ig_ref = {.d = fixedmathftoi(-10.0f / 200.0f), .q = 0, .z = 0},
			.ii_d = {.d = 0, .q = 0, .z = 0},
			.vc_d = {.d = 0, .q = 0, .z = 0},
			.spll_3ph_1 = {.v_q = {0.0f, 0.0f}, .ylf = {0.0f, 0.0f}, .fo = 0.0f, .fn = 50.0f, .theta = {0.0f, 0.0f}, .delta_t = TLVSI_CONFIG_ts, .lpf_coeff = {.b0 = 166.877556f, .b1 = -166.322444f}}
	};

	//float theta, phi;
	float phi;
//    int32_t J, Jk;
//    int32_t co, si;
	fmint_t J, Jk;
	fmint_t co, si;
    uint32_t k;

//    int32_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;
    fmint_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;

    vsi.theta = theta;

    /* Pre-computes sin and cos for DQ0 transforms */
//    si = (int64_t)(TLVSI_CONFIG_FTOI_K * sinf(vsi.theta));
//    co = (int64_t)(TLVSI_CONFIG_FTOI_K * cosf(vsi.theta));

//    tptransformsABCDQ0(ii, &vsi.ii_d, si, co);
//    tptransformsABCDQ0(ig, &vsi.ig_d, si, co);
//    tptransformsABCDQ0(vc, &vsi.vc_d, si, co);
//    tptransformsABCDQ0(vg, &vsi.vg_k, si, co);
//
//    tppllRun(vsi.vg_k.q, &vsi.spll_3ph_1);
//    vsi.theta = vsi.spll_3ph_1.theta[1];

//    vsi.ii_k = *ii;
//    vsi.ig_k = *ig;
//    vsi.vc_k = *vc;
//    vsi.vg_k = *vg;

    vsi.ii_d = *ii;
    vsi.ig_d = *ig;
    vsi.vc_d = *vc;
    vsi.vg_k = *vg;

    /* Delay compensation */
    tlvsiPredictFixed(&vsi.ii_k, &vsi.ii_d, &vsi.ig_k, &vsi.ig_d,
                 &vsi.vc_k, &vsi.vc_d, &vsi.vg_k, vsi.theta, vsi.sw);

    /* References for filter cap. voltage and inverter current */
    vsi.vc_ref.d = fixedmul( TLVSI_CONFIG_k8_int, vsi.ig_ref.q) + vsi.vg_k.d;
    vsi.vc_ref.q = fixedmul(-TLVSI_CONFIG_k8_int, vsi.ig_ref.d) + vsi.vg_k.q;

    vsi.ii_ref.d = vsi.ig_ref.d + fixedmul( TLVSI_CONFIG_k9_int, vsi.vc_ref.q);
    vsi.ii_ref.q = vsi.ig_ref.q + fixedmul(-TLVSI_CONFIG_k9_int, vsi.vc_ref.d);

    /* Predicts for each possible switching combination */
    ii_d_constant = vsi.ii_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, vsi.ii_k.q)
    		+ fixedmul(TLVSI_CONFIG_k2_int, vsi.vc_k.d);

    ii_q_constant = vsi.ii_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, vsi.ii_k.d)
			+ fixedmul(TLVSI_CONFIG_k2_int, vsi.vc_k.q);

    vsi.ii_k_1.d = ii_d_constant;
    vsi.ii_k_1.q = ii_q_constant;

    vc_d_constant = vsi.vc_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, vsi.vc_k.q)
			+ fixedmul(TLVSI_CONFIG_k4_int, vsi.ii_k.d)
			+ fixedmul(TLVSI_CONFIG_k5_int, vsi.ig_k.d);
    vc_q_constant = vsi.vc_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, vsi.vc_k.d)
			+ fixedmul(TLVSI_CONFIG_k4_int, vsi.ii_k.q)
			+ fixedmul(TLVSI_CONFIG_k5_int, vsi.ig_k.q);

    vsi.vc_k_1.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_d_constant);
    vsi.vc_k_1.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_q_constant);

    ig_d_constant = vsi.ig_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, vsi.ig_k.q)
			+ fixedmul(TLVSI_CONFIG_k6_int, vsi.vc_k.d)
			+ fixedmul(TLVSI_CONFIG_k7_int, vsi.vg_k.d);
    ig_q_constant = vsi.ig_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, vsi.ig_k.d)
			+ fixedmul(TLVSI_CONFIG_k6_int, vsi.vc_k.q)
			+ fixedmul(TLVSI_CONFIG_k7_int, vsi.vg_k.q);

    vsi.ig_k_1.d = ig_d_constant + fixedmul(TLVSI_CONFIG_k6_int, vsi.vc_k_1.d);
    vsi.ig_k_1.q = ig_q_constant + fixedmul(TLVSI_CONFIG_k6_int, vsi.vc_k_1.q);

//    tlvsiPredict(vsi, &vsi->ii_k_1, &vsi->ii_k, &vsi->ig_k_1, &vsi->ig_k,
//                 &vsi->vc_k_1, &vsi->vc_k, &vsi.vg_k, vsi->theta, 0);

    Jk = tlvsiCostFixed(&vsi.ii_k_1, &vsi.ii_ref, &vsi.ig_k_1, &vsi.ig_ref,
                   &vsi.vc_k_1, &vsi.vc_ref);
    J = Jk;
    vsi.sw = 0;

    phi = 0;
    for(k = 1; k < 7; k++){

        theta = vsi.theta - phi;
        co = fixedmathftoi( cosf(theta) );
        si = fixedmathftoi( sinf(theta) );

        co = fixedmul(TLVSI_CONFIG_k3_int, co);
        si = fixedmul(TLVSI_CONFIG_k3_int, si);

        vsi.ii_k_1.d = ii_d_constant + co;
        vsi.ii_k_1.q = ii_q_constant - si;

        vsi.vc_k_1.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, vsi.ii_k_1.d);
        vsi.vc_k_1.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, vsi.ii_k_1.q);

        vsi.ig_k_1.d = ig_d_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vsi.vc_k_1.d);
        vsi.ig_k_1.q = ig_q_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vsi.vc_k_1.q);

//        tlvsiPredict(vsi, &vsi->ii_k_1, &vsi->ii_k, &vsi->ig_k_1, &vsi->ig_k,
//                     &vsi->vc_k_1, &vsi->vc_k, &vsi->vg_k, vsi->theta, k);

        Jk = tlvsiCostFixed(&vsi.ii_k_1, &vsi.ii_ref,
                       &vsi.ig_k_1, &vsi.ig_ref, &vsi.vc_k_1, &vsi.vc_ref);

        if(Jk < J){
            J = Jk;
            vsi.sw = k;
        }

        phi += 1.0471975511965976f;
    }

    if( Jopt != 0 ) *Jopt = J;

    return vsi.sw;
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
void tlvsiPredictFixed(psdtypesDQ0int_t *ii_k_1, psdtypesDQ0int_t *ii_k,
				  	   psdtypesDQ0int_t *ig_k_1, psdtypesDQ0int_t *ig_k,
					   psdtypesDQ0int_t *vc_k_1, psdtypesDQ0int_t *vc_k,
					   psdtypesDQ0int_t *vg_k,
					   float theta, uint32_t sw){

	fmint_t co, si;

	ii_k_1->d = ii_k->d + fixedmul(TLVSI_CONFIG_k1_int, ii_k->q) + fixedmul(TLVSI_CONFIG_k2_int, vc_k->d);
	ii_k_1->q = ii_k->q - fixedmul(TLVSI_CONFIG_k1_int, ii_k->d) + fixedmul(TLVSI_CONFIG_k2_int, vc_k->q);
	if(sw != 0){

		theta = theta - (sw-1)*1.0471975511965976f;

		co = fixedmathftoi( cosf(theta) );
        si = fixedmathftoi( sinf(theta) );

        co = fixedmul(TLVSI_CONFIG_k3_int, co);
        si = fixedmul(TLVSI_CONFIG_k3_int, si);

		ii_k_1->d = ii_k_1->d + co;
		ii_k_1->q = ii_k_1->q - si;
	}

	vc_k_1->d = vc_k->d + fixedmul(TLVSI_CONFIG_k1_int, vc_k->q)
			+ fixedmul(TLVSI_CONFIG_k4_int, (ii_k_1->d + ii_k->d) )
			+ fixedmul(TLVSI_CONFIG_k5_int, ig_k->d);
	vc_k_1->q = vc_k->q - fixedmul(TLVSI_CONFIG_k1_int, vc_k->d)
			+ fixedmul(TLVSI_CONFIG_k4_int, (ii_k_1->q + ii_k->q))
			+ fixedmul(TLVSI_CONFIG_k5_int, ig_k->q);

	ig_k_1->d = ig_k->d + fixedmul(TLVSI_CONFIG_k1_int, ig_k->q)
			+ fixedmul(TLVSI_CONFIG_k6_int, (vc_k_1->d + vc_k->d))
			+ fixedmul(TLVSI_CONFIG_k7_int, vg_k->d);
	ig_k_1->q = ig_k->q - fixedmul(TLVSI_CONFIG_k1_int, ig_k->d)
			+ fixedmul(TLVSI_CONFIG_k6_int, (vc_k_1->q + vc_k->q))
			+ fixedmul(TLVSI_CONFIG_k7_int, vg_k->q);
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
fmint_t tlvsiCostFixed(psdtypesDQ0int_t *ii, psdtypesDQ0int_t *ii_ref,
					   psdtypesDQ0int_t *ig, psdtypesDQ0int_t *ig_ref,
					   psdtypesDQ0int_t *vc, psdtypesDQ0int_t *vc_ref){

	fmint_t e_ii[2], e_ig[2], e_vc[2];
	fmint_t J;

	e_ii[0] = ii->d - ii_ref->d;
	e_ii[1] = ii->q - ii_ref->q;

	e_ig[0] = ig->d - ig_ref->d;
	e_ig[1] = ig->q - ig_ref->q;

	e_vc[0] = vc->d - vc_ref->d;
	e_vc[1] = vc->q - vc_ref->q;

	J = fixedmul(TLVSI_CONFIG_w_ii_int, ( fixedmul(e_ii[0], e_ii[0]) + fixedmul(e_ii[1], e_ii[1]) ))
	  + fixedmul(TLVSI_CONFIG_w_ig_int, ( fixedmul(e_ig[0], e_ig[0]) + fixedmul(e_ig[1], e_ig[1]) ))
	  + fixedmul(TLVSI_CONFIG_w_vc_int, ( fixedmul(e_vc[0], e_vc[0]) + fixedmul(e_vc[1], e_vc[1]) ));

	return J;
}
//---------------------------------------------------------------------------
//===========================================================================

