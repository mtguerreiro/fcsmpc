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

#include <stdbool.h>

#include "tpt.h"
//#include <assert.h>

#include "tlvsi_config.h"

int64_t sine(uint64_t value);
int64_t cosine(uint64_t value);
int64_t clamp_overflow(int64_t value, int width);

//===========================================================================

#define TLVSI_CONFIG_k1			((float)(TLVSI_CONFIG_ts * TLVSI_CONFIG_w))
#define TLVSI_CONFIG_k2			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Li))
#define TLVSI_CONFIG_k3			((float)((-TLVSI_CONFIG_ts / TLVSI_CONFIG_Li) * (2.0f/3.0f) * TLVSI_CONFIG_V_dc))
#define TLVSI_CONFIG_k4			((float)(-0.5f * TLVSI_CONFIG_ts / TLVSI_CONFIG_Cf))
#define TLVSI_CONFIG_k5			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Cf))
#define TLVSI_CONFIG_k6			((float)(-0.5f * TLVSI_CONFIG_ts / TLVSI_CONFIG_Lg))
#define TLVSI_CONFIG_k7			((float)(TLVSI_CONFIG_ts / TLVSI_CONFIG_Lg))

#define TLVSI_CONFIG_k8			(TLVSI_CONFIG_w * TLVSI_CONFIG_Lg)
#define TLVSI_CONFIG_k9			(TLVSI_CONFIG_w * TLVSI_CONFIG_Cf)

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

#define TLVSI_II_D_P            TLVSI_CONFIG_II_D_LIM
#define TLVSI_II_D_N           -TLVSI_CONFIG_II_D_LIM

#define TLVSI_II_Q_P            TLVSI_CONFIG_II_Q_LIM
#define TLVSI_II_Q_N           -TLVSI_CONFIG_II_Q_LIM

#define TLVSI_II_D_P_int        fixedmathftoi(TLVSI_CONFIG_II_D_LIM)
#define TLVSI_II_D_N_int        fixedmathftoi(-TLVSI_CONFIG_II_D_LIM)

#define TLVSI_II_Q_P_int        fixedmathftoi(TLVSI_CONFIG_II_Q_LIM)
#define TLVSI_II_Q_N_int        fixedmathftoi(-TLVSI_CONFIG_II_Q_LIM)

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
uint32_t tlvsiOpt(
		tpdtDq0_t *ii, tpdtDq0_t *ig,
		tpdtDq0_t *vc, tpdtDq0_t *vg,
		tpdtDq0_t* ig_ref,
		tpdtfloat_t theta, tpdtfloat_t *Jopt
		){

    static uint32_t sw = 0;
    tpdtfloat_t thetak, phi;
    tpdtfloat_t J, Jk;
    tpdtfloat_t co, si;
    uint32_t k;

    tpdtfloat_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;
    
    tpdtDq0_t ig_k, ig_k_1;
    tpdtDq0_t ii_k, ii_k_1, ii_ref;
    tpdtDq0_t vc_k, vc_k_1, vc_ref;

    /* Delay compensation */
    tlvsiPredict(&ii_k, ii, &ig_k, ig,
                 &vc_k, vc, vg, theta, sw);

    /* References for filter cap. voltage and inverter current */
    vc_ref.d = ( TLVSI_CONFIG_w * TLVSI_CONFIG_Lg) * ig_ref->q + vg->d;
    vc_ref.q = (-TLVSI_CONFIG_w * TLVSI_CONFIG_Lg) * ig_ref->d + vg->q;

    ii_ref.d = ig_ref->d + ( TLVSI_CONFIG_w * TLVSI_CONFIG_Cf) * vc_ref.q;
    ii_ref.q = ig_ref->q + (-TLVSI_CONFIG_w * TLVSI_CONFIG_Cf) * vc_ref.d;

    /* Predicts for each possible switching combination */
    ii_d_constant = ii_k.d + TLVSI_CONFIG_k1 * ii_k.q + TLVSI_CONFIG_k2 * vc_k.d;
    ii_q_constant = ii_k.q - TLVSI_CONFIG_k1 * ii_k.d + TLVSI_CONFIG_k2 * vc_k.q;

    ii_k_1.d = ii_d_constant;
    ii_k_1.q = ii_q_constant;

    vc_d_constant = vc_k.d + TLVSI_CONFIG_k1 * vc_k.q + TLVSI_CONFIG_k4 * ii_k.d + TLVSI_CONFIG_k5 * ig_k.d;
    vc_q_constant = vc_k.q - TLVSI_CONFIG_k1 * vc_k.d + TLVSI_CONFIG_k4 * ii_k.q + TLVSI_CONFIG_k5 * ig_k.q;

    vc_k_1.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_d_constant);
    vc_k_1.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_q_constant);

    ig_d_constant = ig_k.d + TLVSI_CONFIG_k1 * ig_k.q + TLVSI_CONFIG_k6 * vc_k.d + TLVSI_CONFIG_k7 * vg->d;
    ig_q_constant = ig_k.q - TLVSI_CONFIG_k1 * ig_k.d + TLVSI_CONFIG_k6 * vc_k.q + TLVSI_CONFIG_k7 * vg->q;

    ig_k_1.d = ig_d_constant + TLVSI_CONFIG_k6 * (vc_k_1.d);
    ig_k_1.q = ig_q_constant + TLVSI_CONFIG_k6 * (vc_k_1.q);

    if( (ii_k_1.d > TLVSI_II_D_P) || (ii_k_1.d < TLVSI_II_D_N) || (ii_k_1.q > TLVSI_II_Q_P) || (ii_k_1.q < TLVSI_II_Q_N) ){
        Jk = 1e30;
        J = Jk;  
    }
    else{
        Jk = tlvsiCost(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
                       &vc_k_1, &vc_ref);
        J = Jk;  
        sw = 0;
    }
    
    phi = 0;
    for(k = 1; k < 7; k++){

        thetak = theta - phi;
        co = TLVSI_CONFIG_k3 * cosf(thetak);
        si = TLVSI_CONFIG_k3 * sinf(thetak);

        ii_k_1.d = ii_d_constant + co;
        ii_k_1.q = ii_q_constant - si;

        vc_k_1.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_k_1.d);
        vc_k_1.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_k_1.q);

        ig_k_1.d = ig_d_constant + TLVSI_CONFIG_k6 * (vc_k_1.d);
        ig_k_1.q = ig_q_constant + TLVSI_CONFIG_k6 * (vc_k_1.q);

        if( (ii_k_1.d > TLVSI_II_D_P) || (ii_k_1.d < TLVSI_II_D_N) || (ii_k_1.q > TLVSI_II_Q_P) || (ii_k_1.q < TLVSI_II_Q_N) ){
            Jk = 1e30;
        }
        else{
            Jk = tlvsiCost(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
                           &vc_k_1, &vc_ref);    
        }

        if(Jk < J){
            J = Jk;
            sw = k;
        }

        phi += 1.0471975511965976f;
    }

    if( Jopt != 0 ) *Jopt = J;

    return sw;
}
//---------------------------------------------------------------------------
uint32_t tlvsiOpt2(
		tpdtDq0_t *ii, tpdtDq0_t *ig,
		tpdtDq0_t *vc, tpdtDq0_t *vg,
		tpdtDq0_t* ig_ref,
		tpdtfloat_t theta, tpdtfloat_t *Jopt
		){

    static uint32_t sw = 0;
    tpdtfloat_t thetak, phi;
    tpdtfloat_t J, Jk;
    tpdtfloat_t co, si;
    uint32_t k;

    tpdtfloat_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;
    
    tpdtDq0_t ig_k, ig_k_1;
    tpdtDq0_t ii_k, ii_k_1, ii_ref;
    tpdtDq0_t vc_k, vc_k_1, vc_ref;

    /* Delay compensation */
    tlvsiPredict(&ii_k, ii, &ig_k, ig,
                 &vc_k, vc, vg, theta, sw);

    /* References for filter cap. voltage and inverter current */
    vc_ref.d = ( TLVSI_CONFIG_w * TLVSI_CONFIG_Lg) * ig_ref->q + vg->d;
    vc_ref.q = (-TLVSI_CONFIG_w * TLVSI_CONFIG_Lg) * ig_ref->d + vg->q;

    ii_ref.d = ig_ref->d + ( TLVSI_CONFIG_w * TLVSI_CONFIG_Cf) * vc_ref.q;
    ii_ref.q = ig_ref->q + (-TLVSI_CONFIG_w * TLVSI_CONFIG_Cf) * vc_ref.d;

    /* Predicts for each possible switching combination */
//    sw = tlvsiOpt1(&ii_k, &ig_k, &vc_k, vg, &ii_ref, ig_ref, &vc_ref, theta, &J);

    ii_d_constant = ii_k.d + TLVSI_CONFIG_k1 * ii_k.q + TLVSI_CONFIG_k2 * vc_k.d;
    ii_q_constant = ii_k.q - TLVSI_CONFIG_k1 * ii_k.d + TLVSI_CONFIG_k2 * vc_k.q;

    ii_k_1.d = ii_d_constant;
    ii_k_1.q = ii_q_constant;

    vc_d_constant = vc_k.d + TLVSI_CONFIG_k1 * vc_k.q + TLVSI_CONFIG_k4 * ii_k.d + TLVSI_CONFIG_k5 * ig_k.d;
    vc_q_constant = vc_k.q - TLVSI_CONFIG_k1 * vc_k.d + TLVSI_CONFIG_k4 * ii_k.q + TLVSI_CONFIG_k5 * ig_k.q;

    vc_k_1.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_d_constant);
    vc_k_1.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_q_constant);

    ig_d_constant = ig_k.d + TLVSI_CONFIG_k1 * ig_k.q + TLVSI_CONFIG_k6 * vc_k.d + TLVSI_CONFIG_k7 * vg->d;
    ig_q_constant = ig_k.q - TLVSI_CONFIG_k1 * ig_k.d + TLVSI_CONFIG_k6 * vc_k.q + TLVSI_CONFIG_k7 * vg->q;

    ig_k_1.d = ig_d_constant + TLVSI_CONFIG_k6 * (vc_k_1.d);
    ig_k_1.q = ig_q_constant + TLVSI_CONFIG_k6 * (vc_k_1.q);

    //Jk = tlvsiCost(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
    //               &vc_k_1, &vc_ref);
    //J = Jk;
    if( (ii_k_1.d > TLVSI_II_D_P) || (ii_k_1.d < TLVSI_II_D_N) || (ii_k_1.q > TLVSI_II_Q_P) || (ii_k_1.q < TLVSI_II_Q_N) ){
        Jk = 1e30;
        J = Jk;  
    }
    else{
        tlvsiOpt1(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
        J = Jk;
        sw = 0;
    }

    phi = 0;
    for(k = 1; k < 7; k++){

        thetak = theta - phi;
        co = TLVSI_CONFIG_k3 * cosf(thetak);
        si = TLVSI_CONFIG_k3 * sinf(thetak);

        ii_k_1.d = ii_d_constant + co;
        ii_k_1.q = ii_q_constant - si;

        vc_k_1.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_k_1.d);
        vc_k_1.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_k_1.q);

        ig_k_1.d = ig_d_constant + TLVSI_CONFIG_k6 * (vc_k_1.d);
        ig_k_1.q = ig_q_constant + TLVSI_CONFIG_k6 * (vc_k_1.q);

        if( (ii_k_1.d > TLVSI_II_D_P) || (ii_k_1.d < TLVSI_II_D_N) || (ii_k_1.q > TLVSI_II_Q_P) || (ii_k_1.q < TLVSI_II_Q_N) ){
            Jk = 1e30;
        }
        else{
            tlvsiOpt1(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
        }

//        tlvsiOpt1(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
        //Jk = tlvsiCost(&ii_k_1, &ii_ref,
        //               &ig_k_1, ig_ref, &vc_k_1, &vc_ref);

        if(Jk < J){
            J = Jk;
            sw = k;
        }

        phi += 1.0471975511965976f;
    }

    if( Jopt != 0 ) *Jopt = J;

    return sw;
}
//---------------------------------------------------------------------------
uint32_t tlvsiOpt1(
		tpdtDq0_t *ii, tpdtDq0_t *ig,
		tpdtDq0_t *vc, tpdtDq0_t *vg,
		tpdtDq0_t *ii_ref, tpdtDq0_t *ig_ref,
		tpdtDq0_t *vc_ref,
		tpdtfloat_t theta, tpdtfloat_t *Jopt
		){

    static uint32_t sw = 0;
    tpdtfloat_t thetak, phi;
    tpdtfloat_t J, Jk;
    tpdtfloat_t co, si;
    uint32_t k;

    tpdtfloat_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;
    
    tpdtDq0_t ig_k;
    tpdtDq0_t ii_k;
    tpdtDq0_t vc_k;

    /* Predicts for each possible switching combination */
    ii_d_constant = ii->d + TLVSI_CONFIG_k1 * ii->q + TLVSI_CONFIG_k2 * vc->d;
    ii_q_constant = ii->q - TLVSI_CONFIG_k1 * ii->d + TLVSI_CONFIG_k2 * vc->q;

    ii_k.d = ii_d_constant;
    ii_k.q = ii_q_constant;

    vc_d_constant = vc->d + TLVSI_CONFIG_k1 * vc->q + TLVSI_CONFIG_k4 * ii->d + TLVSI_CONFIG_k5 * ig->d;
    vc_q_constant = vc->q - TLVSI_CONFIG_k1 * vc->d + TLVSI_CONFIG_k4 * ii->q + TLVSI_CONFIG_k5 * ig->q;

    vc_k.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_d_constant);
    vc_k.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_q_constant);

    ig_d_constant = ig->d + TLVSI_CONFIG_k1 * ig->q + TLVSI_CONFIG_k6 * vc->d + TLVSI_CONFIG_k7 * vg->d;
    ig_q_constant = ig->q - TLVSI_CONFIG_k1 * ig->d + TLVSI_CONFIG_k6 * vc->q + TLVSI_CONFIG_k7 * vg->q;

    ig_k.d = ig_d_constant + TLVSI_CONFIG_k6 * (vc_k.d);
    ig_k.q = ig_q_constant + TLVSI_CONFIG_k6 * (vc_k.q);

/*
    if( (ii_k.d > TLVSI_II_D_P) || (ii_k.d < TLVSI_II_D_N) || (ii_k.q > TLVSI_II_Q_P) || (ii_k.q < TLVSI_II_Q_N) ){
        Jk = 1e30;
        J = Jk;  
    }
    else{
        Jk = tlvsiCost(&ii_k, ii_ref, &ig_k, ig_ref,
                       &vc_k, vc_ref);
        J = Jk;  
        sw = 0;
    }
*/
    Jk = tlvsiCost(&ii_k, ii_ref, &ig_k, ig_ref,
                   &vc_k, vc_ref);
    J = Jk;
    sw = 0;

    phi = 0;
    for(k = 1; k < 7; k++){

        thetak = theta - phi;
        co = TLVSI_CONFIG_k3 * cosf(thetak);
        si = TLVSI_CONFIG_k3 * sinf(thetak);

        ii_k.d = ii_d_constant + co;
        ii_k.q = ii_q_constant - si;

        vc_k.d = vc_d_constant + TLVSI_CONFIG_k4 * (ii_k.d);
        vc_k.q = vc_q_constant + TLVSI_CONFIG_k4 * (ii_k.q);

        ig_k.d = ig_d_constant + TLVSI_CONFIG_k6 * (vc_k.d);
        ig_k.q = ig_q_constant + TLVSI_CONFIG_k6 * (vc_k.q);
/*
        if( (ii_k.d > TLVSI_II_D_P) || (ii_k.d < TLVSI_II_D_N) || (ii_k.q > TLVSI_II_Q_P) || (ii_k.q < TLVSI_II_Q_N) ){
            Jk = 1e30;
        }
        else{
            Jk = tlvsiCost(&ii_k, ii_ref, &ig_k, ig_ref,
                           &vc_k, vc_ref);
        }
*/   
        Jk = tlvsiCost(&ii_k, ii_ref,
                       &ig_k, ig_ref, &vc_k, vc_ref);

        if(Jk < J){
            J = Jk;
            sw = k;
        }

        phi += 1.0471975511965976f;
    }

    if( Jopt != 0 ) *Jopt = J;

    return sw;
}
//---------------------------------------------------------------------------
uint32_t tlvsiOptFixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ig,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vg,
		tpdtDq0FixedPoint_t * ig_ref,
		tpdtint_t theta, tpdtint_t *Jopt
		){

    static uint32_t sw = 0;
    tpdtint_t thetak, phi;
    tpdtint_t co, si;
    tpdtint_t J, Jk;
    uint32_t k;

    tpdtDq0FixedPoint_t ig_k, ig_k_1;
    tpdtDq0FixedPoint_t ii_k, ii_k_1, ii_ref;
    tpdtDq0FixedPoint_t vc_k, vc_k_1, vc_ref;
    
    tpdtint_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;

    /* Delay compensation */
    tlvsiPredictFixed(&ii_k, ii, &ig_k, ig,
                 &vc_k, vc, vg, theta, sw);

    /* References for filter cap. voltage and inverter current */
    vc_ref.d = fixedmul( TLVSI_CONFIG_k8_int, ig_ref->q) + vg->d;
    vc_ref.q = fixedmul(-TLVSI_CONFIG_k8_int, ig_ref->d) + vg->q;

    ii_ref.d = ig_ref->d + fixedmul( TLVSI_CONFIG_k9_int, vc_ref.q);
    ii_ref.q = ig_ref->q + fixedmul(-TLVSI_CONFIG_k9_int, vc_ref.d);

    /* Predicts for each possible switching combination */
    ii_d_constant = ii_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, ii_k.q)
    		+ fixedmul(TLVSI_CONFIG_k2_int, vc_k.d);

    ii_q_constant = ii_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, ii_k.d)
			+ fixedmul(TLVSI_CONFIG_k2_int, vc_k.q);

    ii_k_1.d = ii_d_constant;
    ii_k_1.q = ii_q_constant;

    vc_d_constant = vc_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, vc_k.q)
			+ fixedmul(TLVSI_CONFIG_k4_int, ii_k.d)
			+ fixedmul(TLVSI_CONFIG_k5_int, ig_k.d);
    vc_q_constant = vc_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, vc_k.d)
			+ fixedmul(TLVSI_CONFIG_k4_int, ii_k.q)
			+ fixedmul(TLVSI_CONFIG_k5_int, ig_k.q);

    vc_k_1.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_d_constant);
    vc_k_1.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_q_constant);

    ig_d_constant = ig_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, ig_k.q)
			+ fixedmul(TLVSI_CONFIG_k6_int, vc_k.d)
			+ fixedmul(TLVSI_CONFIG_k7_int, vg->d);
    ig_q_constant = ig_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, ig_k.d)
			+ fixedmul(TLVSI_CONFIG_k6_int, vc_k.q)
			+ fixedmul(TLVSI_CONFIG_k7_int, vg->q);

    ig_k_1.d = ig_d_constant + fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.d);
    ig_k_1.q = ig_q_constant + fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.q);


    if( (ii_k_1.d > TLVSI_II_D_P_int) || (ii_k_1.d < TLVSI_II_D_N_int) || (ii_k_1.q > TLVSI_II_Q_P_int) || (ii_k_1.q < TLVSI_II_Q_N_int) ){
        Jk = (0xFFFFFFFF >> 1);
        J = Jk;  
    }
    else{
        Jk = tlvsiCostFixed(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
                       &vc_k_1, &vc_ref);
        J = Jk;  
        sw = 0;
    }
  
/*
    Jk = tlvsiCostFixed(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
                   &vc_k_1, &vc_ref);
    J = Jk;
    sw = 0;
*/

    phi = 0;
    for(k = 1; k < 7; k++){

        thetak = theta - phi;

        co = cosine((uint64_t)(fixedmul(thetak, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        co = co << (FIXED_MATH_Q - 18 + 1);

        si = sine((uint64_t)(fixedmul(thetak, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        si = si << (FIXED_MATH_Q - 18 + 1);

        co = fixedmul(TLVSI_CONFIG_k3_int, co);
        si = fixedmul(TLVSI_CONFIG_k3_int, si);

        ii_k_1.d = ii_d_constant + co;
        ii_k_1.q = ii_q_constant - si;

        vc_k_1.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_k_1.d);
        vc_k_1.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_k_1.q);

        ig_k_1.d = ig_d_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.d);
        ig_k_1.q = ig_q_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.q);


        if( (ii_k_1.d > TLVSI_II_D_P_int) || (ii_k_1.d < TLVSI_II_D_N_int) || (ii_k_1.q > TLVSI_II_Q_P_int) || (ii_k_1.q < TLVSI_II_Q_N_int) ){
            Jk = (0xFFFFFFFF >> 1);
        }
        else{
            Jk = tlvsiCostFixed(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
                           &vc_k_1, &vc_ref);
        }

/*
        Jk = tlvsiCostFixed(&ii_k_1, &ii_ref,
                       &ig_k_1, ig_ref, &vc_k_1, &vc_ref);
*/

        if(Jk < J){
            J = Jk;
            sw = k;
        }

        phi += fixedmathftoi(1.0471975511965976f);
    }

    if( Jopt != 0 ) *Jopt = J;

    return sw;
}
//---------------------------------------------------------------------------
uint32_t tlvsiOpt2Fixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ig,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vg,
		tpdtDq0FixedPoint_t* ig_ref,
		tpdtint_t theta, tpdtint_t *Jopt
		){

    static uint32_t sw = 0;
    tpdtint_t thetak, phi;
    tpdtint_t co, si;
    tpdtint_t J, Jk;
    uint32_t k;

    tpdtDq0FixedPoint_t ig_k, ig_k_1;
    tpdtDq0FixedPoint_t ii_k, ii_k_1, ii_ref;
    tpdtDq0FixedPoint_t vc_k, vc_k_1, vc_ref;
    
    tpdtint_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;

    /* Delay compensation */
    tlvsiPredictFixed(&ii_k, ii, &ig_k, ig,
                 &vc_k, vc, vg, theta, sw);

    /* References for filter cap. voltage and inverter current */
    vc_ref.d = fixedmul( TLVSI_CONFIG_k8_int, ig_ref->q) + vg->d;
    vc_ref.q = fixedmul(-TLVSI_CONFIG_k8_int, ig_ref->d) + vg->q;

    ii_ref.d = ig_ref->d + fixedmul( TLVSI_CONFIG_k9_int, vc_ref.q);
    ii_ref.q = ig_ref->q + fixedmul(-TLVSI_CONFIG_k9_int, vc_ref.d);

    /* Predicts for each possible switching combination */
    //sw = tlvsiOpt1Fixed(&ii_k, &ig_k, &vc_k, vg, &ii_ref, ig_ref, &vc_ref, theta, &J);
    
    ii_d_constant = ii_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, ii_k.q)
    		+ fixedmul(TLVSI_CONFIG_k2_int, vc_k.d);

    ii_q_constant = ii_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, ii_k.d)
			+ fixedmul(TLVSI_CONFIG_k2_int, vc_k.q);

    ii_k_1.d = ii_d_constant;
    ii_k_1.q = ii_q_constant;

    vc_d_constant = vc_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, vc_k.q)
			+ fixedmul(TLVSI_CONFIG_k4_int, ii_k.d)
			+ fixedmul(TLVSI_CONFIG_k5_int, ig_k.d);
    vc_q_constant = vc_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, vc_k.d)
			+ fixedmul(TLVSI_CONFIG_k4_int, ii_k.q)
			+ fixedmul(TLVSI_CONFIG_k5_int, ig_k.q);

    vc_k_1.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_d_constant);
    vc_k_1.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_q_constant);

    ig_d_constant = ig_k.d
    		+ fixedmul(TLVSI_CONFIG_k1_int, ig_k.q)
			+ fixedmul(TLVSI_CONFIG_k6_int, vc_k.d)
			+ fixedmul(TLVSI_CONFIG_k7_int, vg->d);
    ig_q_constant = ig_k.q
    		- fixedmul(TLVSI_CONFIG_k1_int, ig_k.d)
			+ fixedmul(TLVSI_CONFIG_k6_int, vc_k.q)
			+ fixedmul(TLVSI_CONFIG_k7_int, vg->q);

    ig_k_1.d = ig_d_constant + fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.d);
    ig_k_1.q = ig_q_constant + fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.q);

    //Jk = tlvsiCostFixed(&ii_k_1, &ii_ref, &ig_k_1, ig_ref,
    //               &vc_k_1, &vc_ref);
    
    if( (ii_k_1.d > TLVSI_II_D_P_int) || (ii_k_1.d < TLVSI_II_D_N_int) || (ii_k_1.q > TLVSI_II_Q_P_int) || (ii_k_1.q < TLVSI_II_Q_N_int) ){
        Jk = (0xFFFFFFFF >> 1);
        J = Jk;  
    }
    else{
        tlvsiOpt1Fixed(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
        J = Jk;  
        sw = 0;
    }
    
//    tlvsiOpt1Fixed(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
    
    J = Jk;
    sw = 0;

    phi = 0;
    for(k = 1; k < 7; k++){

        thetak = theta - phi;

        co = cosine((uint64_t)(fixedmul(thetak, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        co = co << (FIXED_MATH_Q - 18 + 1);

        si = sine((uint64_t)(fixedmul(thetak, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        si = si << (FIXED_MATH_Q - 18 + 1);

        co = fixedmul(TLVSI_CONFIG_k3_int, co);
        si = fixedmul(TLVSI_CONFIG_k3_int, si);

        ii_k_1.d = ii_d_constant + co;
        ii_k_1.q = ii_q_constant - si;

        vc_k_1.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_k_1.d);
        vc_k_1.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_k_1.q);

        ig_k_1.d = ig_d_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.d);
        ig_k_1.q = ig_q_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vc_k_1.q);

        if( (ii_k_1.d > TLVSI_II_D_P_int) || (ii_k_1.d < TLVSI_II_D_N_int) || (ii_k_1.q > TLVSI_II_Q_P_int) || (ii_k_1.q < TLVSI_II_Q_N_int) ){
            Jk = (0xFFFFFFFF >> 1);
        }
        else{
            tlvsiOpt1Fixed(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
        }
    
//        tlvsiOpt1Fixed(&ii_k_1, &ig_k_1, &vc_k_1, vg, &ii_ref, ig_ref, &vc_ref, theta, &Jk);
    
        //Jk = tlvsiCostFixed(&ii_k_1, &ii_ref,
        //               &ig_k_1, ig_ref, &vc_k_1, &vc_ref);

        if(Jk < J){
            J = Jk;
            sw = k;
        }

        phi += fixedmathftoi(1.0471975511965976f);
    }

    if( Jopt != 0 ) *Jopt = J;

    return sw;
}
//---------------------------------------------------------------------------
uint32_t tlvsiOpt1Fixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ig,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vg,
		tpdtDq0FixedPoint_t* ii_ref, tpdtDq0FixedPoint_t* ig_ref,
		tpdtDq0FixedPoint_t* vc_ref,
		tpdtint_t theta, tpdtint_t *Jopt
		){

    static uint32_t sw = 0;
    tpdtint_t thetak, phi;
    tpdtint_t co, si;
    tpdtint_t J, Jk;
    uint32_t k;

    tpdtDq0FixedPoint_t ig_k;
    tpdtDq0FixedPoint_t ii_k;
    tpdtDq0FixedPoint_t vc_k;
    
    tpdtint_t ii_d_constant, ii_q_constant, ig_d_constant, ig_q_constant, vc_d_constant, vc_q_constant;

    ii_d_constant = ii->d
    		+ fixedmul(TLVSI_CONFIG_k1_int, ii->q)
    		+ fixedmul(TLVSI_CONFIG_k2_int, vc->d);

    ii_q_constant = ii->q
    		- fixedmul(TLVSI_CONFIG_k1_int, ii->d)
			+ fixedmul(TLVSI_CONFIG_k2_int, vc->q);

    ii_k.d = ii_d_constant;
    ii_k.q = ii_q_constant;

    vc_d_constant = vc->d
    		+ fixedmul(TLVSI_CONFIG_k1_int, vc->q)
			+ fixedmul(TLVSI_CONFIG_k4_int, ii->d)
			+ fixedmul(TLVSI_CONFIG_k5_int, ig->d);
    vc_q_constant = vc->q
    		- fixedmul(TLVSI_CONFIG_k1_int, vc->d)
			+ fixedmul(TLVSI_CONFIG_k4_int, ii->q)
			+ fixedmul(TLVSI_CONFIG_k5_int, ig->q);

    vc_k.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_d_constant);
    vc_k.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_q_constant);

    ig_d_constant = ig->d
    		+ fixedmul(TLVSI_CONFIG_k1_int, ig->q)
			+ fixedmul(TLVSI_CONFIG_k6_int, vc->d)
			+ fixedmul(TLVSI_CONFIG_k7_int, vg->d);
    ig_q_constant = ig->q
    		- fixedmul(TLVSI_CONFIG_k1_int, ig->d)
			+ fixedmul(TLVSI_CONFIG_k6_int, vc->q)
			+ fixedmul(TLVSI_CONFIG_k7_int, vg->q);

    ig_k.d = ig_d_constant + fixedmul(TLVSI_CONFIG_k6_int, vc_k.d);
    ig_k.q = ig_q_constant + fixedmul(TLVSI_CONFIG_k6_int, vc_k.q);


    if( (ii_k.d > TLVSI_II_D_P_int) || (ii_k.d < TLVSI_II_D_N_int) || (ii_k.q > TLVSI_II_Q_P_int) || (ii_k.q < TLVSI_II_Q_N_int) ){
        Jk = (0xFFFFFFFF >> 1);
        J = Jk;  
    }
    else{
        Jk = tlvsiCostFixed(&ii_k, ii_ref, &ig_k, ig_ref,
                       &vc_k, vc_ref);
        J = Jk;  
        sw = 0;
    }

/*
    Jk = tlvsiCostFixed(&ii_k, ii_ref, &ig_k, ig_ref,
                   &vc_k, vc_ref);
    J = Jk;
    sw = 0;
*/

    phi = 0;
    for(k = 1; k < 7; k++){

        thetak = theta - phi;

        co = cosine((uint64_t)(fixedmul(thetak, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        co = co << (FIXED_MATH_Q - 18 + 1);

        si = sine((uint64_t)(fixedmul(thetak, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        si = si << (FIXED_MATH_Q - 18 + 1);

        co = fixedmul(TLVSI_CONFIG_k3_int, co);
        si = fixedmul(TLVSI_CONFIG_k3_int, si);

        ii_k.d = ii_d_constant + co;
        ii_k.q = ii_q_constant - si;

        vc_k.d = vc_d_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_k.d);
        vc_k.q = vc_q_constant + fixedmul(TLVSI_CONFIG_k4_int, ii_k.q);

        ig_k.d = ig_d_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vc_k.d);
        ig_k.q = ig_q_constant
        		+ fixedmul(TLVSI_CONFIG_k6_int, vc_k.q);

        if( (ii_k.d > TLVSI_II_D_P_int) || (ii_k.d < TLVSI_II_D_N_int) || (ii_k.q > TLVSI_II_Q_P_int) || (ii_k.q < TLVSI_II_Q_N_int) ){
            Jk = (0xFFFFFFFF >> 1);
        }
        else{
            Jk = tlvsiCostFixed(&ii_k, ii_ref, &ig_k, ig_ref,
                           &vc_k, vc_ref);
        }

/*
        Jk = tlvsiCostFixed(&ii_k, ii_ref,
                       &ig_k, ig_ref, &vc_k, vc_ref);
*/

        if(Jk < J){
            J = Jk;
            sw = k;
        }

        phi += fixedmathftoi(1.0471975511965976f);
    }

    if( Jopt != 0 ) *Jopt = J;
    
    return sw;
}
//---------------------------------------------------------------------------
void tlvsiPredict(
		tpdtDq0_t *ii_k_1, tpdtDq0_t *ii_k,
		tpdtDq0_t *ig_k_1, tpdtDq0_t *ig_k,
		tpdtDq0_t *vc_k_1, tpdtDq0_t *vc_k,
		tpdtDq0_t *vg_k,
		tpdtfloat_t theta, uint32_t sw
		){

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
void tlvsiPredictFixed(
		tpdtDq0FixedPoint_t *ii_k_1, tpdtDq0FixedPoint_t *ii_k,
		tpdtDq0FixedPoint_t *ig_k_1, tpdtDq0FixedPoint_t *ig_k,
		tpdtDq0FixedPoint_t *vc_k_1, tpdtDq0FixedPoint_t *vc_k,
		tpdtDq0FixedPoint_t *vg_k,
		tpdtint_t theta, uint32_t sw
		){

	tpdtint_t co, si;

	ii_k_1->d = ii_k->d + fixedmul(TLVSI_CONFIG_k1_int, ii_k->q) + fixedmul(TLVSI_CONFIG_k2_int, vc_k->d);
	ii_k_1->q = ii_k->q - fixedmul(TLVSI_CONFIG_k1_int, ii_k->d) + fixedmul(TLVSI_CONFIG_k2_int, vc_k->q);
	if(sw != 0){

		theta = theta - fixedmul( ( (sw-1) << FIXED_MATH_Q ), ( fixedmathftoi(1.0471975511965976f) ) );

        co = cosine((uint64_t)(fixedmul(theta, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        co = co << (FIXED_MATH_Q - 18 + 1);

        si = sine((uint64_t)(fixedmul(theta, (fixedmathftoi(0.15915494309189535f))) >> (FIXED_MATH_Q - 20)));
        si = si << (FIXED_MATH_Q - 18 + 1);

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
tpdtfloat_t tlvsiCost(
		tpdtDq0_t *ii, tpdtDq0_t *ii_ref,
		tpdtDq0_t *ig, tpdtDq0_t *ig_ref,
		tpdtDq0_t *vc, tpdtDq0_t *vc_ref
		){

	tpdtfloat_t e_ii[2], e_ig[2], e_vc[2];
	tpdtfloat_t J;

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
tpdtint_t tlvsiCostFixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ii_ref,
		tpdtDq0FixedPoint_t *ig, tpdtDq0FixedPoint_t *ig_ref,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vc_ref
		){

	tpdtint_t e_ii[2], e_ig[2], e_vc[2];
	tpdtint_t J;

	e_ii[0] = ii->d - ii_ref->d;
	e_ii[1] = ii->q - ii_ref->q;

	e_ig[0] = ig->d - ig_ref->d;
	e_ig[1] = ig->q - ig_ref->q;

	e_vc[0] = vc->d - vc_ref->d;
	e_vc[1] = vc->q - vc_ref->q;

	J = fixedmul(TLVSI_CONFIG_w_ii_int,  fixedmul(e_ii[0], e_ii[0]) + fixedmul(e_ii[1], e_ii[1]) )
	  + fixedmul(TLVSI_CONFIG_w_ig_int,  fixedmul(e_ig[0], e_ig[0]) + fixedmul(e_ig[1], e_ig[1]) )
	  + fixedmul(TLVSI_CONFIG_w_vc_int,  fixedmul(e_vc[0], e_vc[0]) + fixedmul(e_vc[1], e_vc[1]) );

	return J;
}
//---------------------------------------------------------------------------
//===========================================================================

/*
 * Implementation of the fixed-point sine and cosine functions were extracted
 * from: https://github.com/dsnet/tri-approx
 *
 * There were only very small changes.
 */

// Fixed-point sine approximation. Normalized for an input domain of [0,1)
// instead of the usual domain of [0,2*PI).
//
// Uses Taylor series approximation for sine centered at zero:
//  sine(2*PI*x) = 0 + (2*PI*x)^1/1! - (2*PI*x)^3/3!
//                   + (2*PI*x)^5/5! - (2*PI*x)^7/7!
//               = k_1*x^1 - k_3*x^3 + k_5*x^5 - k_7*x^7
//
// The bit-width of 18 appears often because it is the width of hardware
// multipliers on Altera FPGAs.
//
// Input: 20-bit unsigned fixed point integer upscaled by 2^20
// Output: 18-bit two's complement fixed point integer upscaled by 2^17
int64_t sine(uint64_t value) {
    // These are polynomial constants generated for each term in the Taylor
    // series. They have been upscaled to the largest value that fits within
    // 18-bits for greatest precision. The constants labeled with [ADJ] have
    // been manually adjusted to increase accuracy.
    const uint64_t k1 = 205887; // k1 = round((2*PI)^1/1! * 2^15)
    const uint64_t k3 = 169336; // k3 = round((2*PI)^3/3! * 2^12)
    const uint64_t k5 = 167014; // k5 = round((2*PI)^5/5! * 2^11) [ADJ]
    const uint64_t k7 = 150000; // k7 = round((2*PI)^7/7! * 2^11) [ADJ]

    // Uses symmetric properties of sine to get more accurate results
    // Normalize the x value to a 18-bit value upscaled by 2^20
    bool high0 = ((value >> 19) & 0x01);
    bool high1 = ((value >> 18) & 0x01);
    uint64_t x1 = value & 0x3ffff; // Strip two highest bits
    if (high1) {
        x1 = (((uint64_t)0x1 << 18) - x1) & 0x3ffff;
    }
    bool negative = high0;
    bool one = (x1 == 0) && high1;

    // Compute the power values (most of these must be done in series)
    uint64_t x2 = ((x1 * x1) >> 18); // Scale: 2^22
    uint64_t x3 = ((x2 * x1) >> 18); // Scale: 2^24
    uint64_t x5 = ((x2 * x3) >> 18); // Scale: 2^28
    uint64_t x7 = ((x2 * x5) >> 18); // Scale: 2^32

    // Compute the polynomial values (these can be done in parallel)
    uint64_t kx1 = ((k1 * x1) >> 17); // Scale: 2^18
    uint64_t kx3 = ((k3 * x3) >> 18); // Scale: 2^18
    uint64_t kx5 = ((k5 * x5) >> 21); // Scale: 2^18
    uint64_t kx7 = ((k7 * x7) >> 25); // Scale: 2^18

    // Add all the terms together (these can be done in series-parallel)
    int64_t sum = kx1 - kx3 + kx5 - kx7; // Scale: 2^18
    sum = sum >> 1; // Scale: 2^17

    // Perform reflection math and corrections
    if (one) { // Check if sum should be one
        sum = (1 << 17);
    }
    if (negative) { // Check if the sum should be negative
        sum = ~sum + 1;
    }
    return clamp_overflow(sum, 18);
}

// Fixed-point cosine approximation. Normalized for an input domain of [0,1)
// instead of the usual domain of [0,2*PI).
//
// Uses Taylor series approximation for cosines centered at zero:
//  cosine(2*PI*x) = 1 - (2*PI*x)^2/2! + (2*PI*x)^4/4!
//                     - (2*PI*x)^6/6! + (2*PI*x)^8/8!
//                 = 1 - k_2*x^2 + k_4*x^4 - k_6*x^6 + k_8*x^8
//
// The bit-width of 18 appears often because it is the width of hardware
// multipliers on Altera FPGAs.
//
// Input: 20-bit unsigned fixed point integer upscaled by 2^20
// Output: 18-bit two's complement fixed point integer upscaled by 2^17
int64_t cosine(uint64_t value) {
    // These are polynomial constants generated for each term in the Taylor
    // series. They have been upscaled to the largest value that fits within
    // 18-bits for greatest precision. The constants labeled with [ADJ] have
    // been manually adjusted to increase accuracy.
    const uint64_t k2 = 161704; // k2 = round((2*PI)^2/2! * 2^13)
    const uint64_t k4 = 132996; // k4 = round((2*PI)^4/4! * 2^11)
    const uint64_t k6 = 175016; // k6 = round((2*PI)^6/6! * 2^11)
    const uint64_t k8 = 241700; // k8 = round((2*PI)^8/8! * 2^12) [ADJ]

    // Uses symmetric properties of cosine to get more accurate results
    // Normalize the x value to a 18-bit value upscaled by 2^20
    bool high0 = ((value >> 19) & 0x01);
    bool high1 = ((value >> 18) & 0x01);
    uint64_t x1 = value & 0x3ffff; // Strip two highest bits
    if (high1) {
        x1 = (((uint64_t)0x1 << 18) - x1) & 0x3ffff;
    }
    bool negative = high0 ^ high1;
    bool zero = (x1 == 0) && high1;

    // Compute the power values (most of these must be done in series)
    uint64_t x2 = ((x1 * x1) >> 18); // Scale: 2^22
    uint64_t x4 = ((x2 * x2) >> 18); // Scale: 2^26
    uint64_t x6 = ((x4 * x2) >> 18); // Scale: 2^30
    uint64_t x8 = ((x4 * x4) >> 18); // Scale: 2^34

    // Compute the polynomial values (these can be done in parallel)
    uint64_t kx2 = ((k2 * x2) >> 17); // Scale: 2^18
    uint64_t kx4 = ((k4 * x4) >> 19); // Scale: 2^18
    uint64_t kx6 = ((k6 * x6) >> 23); // Scale: 2^18
    uint64_t kx8 = ((k8 * x8) >> 28); // Scale: 2^18

    // Add all the terms together (these can be done in series-parallel)
    int64_t sum = ((int64_t)1 << 18) - kx2 + kx4 - kx6 + kx8; // Scale: 2^18
    sum = sum >> 1; // Scale: 2^17

    // Perform reflection math and corrections
    if (zero) { // Check if sum should be zero
        sum = 0;
    }
    if (negative) { // Check if the sum should be negative
        sum = ~sum + 1;
    }
    return clamp_overflow(sum, 18);
}

// Clamp an overflowed fixed-point two's complement value.
int64_t clamp_overflow(int64_t value, int width) {
    bool high0 = ((value >> (width+0)) & 0x01);
    bool high1 = ((value >> (width-1)) & 0x01);
    if (high0 ^ high1) {
        if (high0) {
            value = (-1 << (width-1));
        } else {
            value = (1 << (width-1)) - 1;
        }
    }

//    // Sanity check to ensure there is no overflow
//    int64_t high = value >> (width-1);
//    assert((high == 0x0000000000000000LL) || (high == 0xFFFFFFFFFFFFFFFFLL));
    return value;
}
