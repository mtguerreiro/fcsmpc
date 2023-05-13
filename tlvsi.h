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

#include "psutils/pll.h"
#include "psutils/tpdt.h"

#include "psutils/fixedmath.h"
//===========================================================================

//===========================================================================
/*------------------------------ Definitions ------------------------------*/
//===========================================================================

typedef struct{

	tpdtAbc_t ig_abc;
	tpdtDq0_t ig_d;
	tpdtDq0_t ig_k;
	tpdtDq0_t ig_k_1;
	tpdtDq0_t ig_ref;

	tpdtAbc_t ii_abc;
    tpdtDq0_t ii_d;
    tpdtDq0_t ii_k;
    tpdtDq0_t ii_k_1;
    tpdtDq0_t ii_ref;

    tpdtAbc_t vc_abc;
    tpdtDq0_t vc_d;
    tpdtDq0_t vc_k;
    tpdtDq0_t vc_k_1;
    tpdtDq0_t vc_ref;

    tpdtAbc_t vg_abc;
    tpdtDq0_t vg_k;

    float theta;

    uint32_t sw;

    pll_t spll_3ph_1;

} tlvsiLCLPredictData_t;

typedef struct{

	tpdtAbc_t ig_abc;
	tpdtDq0FixedPoint_t ig_d;
	tpdtDq0FixedPoint_t ig_k;
	tpdtDq0FixedPoint_t ig_k_1;
	tpdtDq0FixedPoint_t ig_ref;

    tpdtAbc_t ii_abc;
    tpdtDq0FixedPoint_t ii_d;
    tpdtDq0FixedPoint_t ii_k;
    tpdtDq0FixedPoint_t ii_k_1;
    tpdtDq0FixedPoint_t ii_ref;

    tpdtAbc_t vc_abc;
    tpdtDq0FixedPoint_t vc_d;
    tpdtDq0FixedPoint_t vc_k;
    tpdtDq0FixedPoint_t vc_k_1;
    tpdtDq0FixedPoint_t vc_ref;

    tpdtAbc_t vg_abc;
    tpdtDq0FixedPoint_t vg_k;

    fmint_t theta;

    uint32_t sw;

    pllFixedPoint_t spll_3ph_1;

} tlvsiLCLPredictDataInt_t;
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
uint32_t tlvsiOpt(
		tpdtDq0_t *ii, tpdtDq0_t *ig,
		tpdtDq0_t *vc, tpdtDq0_t *vg,
		tpdtDq0_t* ig_ref,
		tpdtfloat_t theta, tpdtfloat_t *Jopt
		);

void tlvsiPredict(
		tpdtDq0_t *ii_k_1, tpdtDq0_t *ii_k,
		tpdtDq0_t *ig_k_1, tpdtDq0_t *ig_k,
		tpdtDq0_t *vc_k_1, tpdtDq0_t *vc_k,
		tpdtDq0_t *vg_k,
		tpdtfloat_t theta, uint32_t sw
		);

void tlvsiPredictFixed(
		tpdtDq0FixedPoint_t *ii_k_1, tpdtDq0FixedPoint_t *ii_k,
		tpdtDq0FixedPoint_t *ig_k_1, tpdtDq0FixedPoint_t *ig_k,
		tpdtDq0FixedPoint_t *vc_k_1, tpdtDq0FixedPoint_t *vc_k,
		tpdtDq0FixedPoint_t *vg_k,
		tpdtint_t theta, uint32_t sw
		);

uint32_t tlvsiOptFixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ig,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vg,
		tpdtDq0FixedPoint_t * ig_ref,
		tpdtint_t theta, tpdtint_t *Jopt
		);

tpdtfloat_t tlvsiCost(
		tpdtDq0_t *ii, tpdtDq0_t *ii_ref,
		tpdtDq0_t *ig, tpdtDq0_t *ig_ref,
		tpdtDq0_t *vc, tpdtDq0_t *vc_ref
		);

tpdtint_t tlvsiCostFixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ii_ref,
		tpdtDq0FixedPoint_t *ig, tpdtDq0FixedPoint_t *ig_ref,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vc_ref
		);
                     
uint32_t tlvsiOpt2Fixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ig,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vg,
		tpdtDq0FixedPoint_t* ig_ref,
		tpdtint_t theta, tpdtint_t *Jopt
		);

uint32_t tlvsiOpt1Fixed(
		tpdtDq0FixedPoint_t *ii, tpdtDq0FixedPoint_t *ig,
		tpdtDq0FixedPoint_t *vc, tpdtDq0FixedPoint_t *vg,
		tpdtDq0FixedPoint_t* ii_ref, tpdtDq0FixedPoint_t* ig_ref,
		tpdtDq0FixedPoint_t* vc_ref,
		tpdtint_t theta, tpdtint_t *Jopt
		);

uint32_t tlvsiOpt2(
		tpdtDq0_t *ii, tpdtDq0_t *ig,
		tpdtDq0_t *vc, tpdtDq0_t *vg,
		tpdtDq0_t* ig_ref,
		tpdtfloat_t theta, tpdtfloat_t *Jopt
		);

uint32_t tlvsiOpt1(
		tpdtDq0_t *ii, tpdtDq0_t *ig,
		tpdtDq0_t *vc, tpdtDq0_t *vg,
		tpdtDq0_t *ii_ref, tpdtDq0_t *ig_ref,
		tpdtDq0_t *vc_ref,
		tpdtfloat_t theta, tpdtfloat_t *Jopt
		);
//===========================================================================



#endif /* TLVSI_H_ */
