/*
 * pll.c
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco.
 *
 */

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include "pll.h"

#include "fixedmath.h"
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
void pllInit(pll_t *pll, pllfloat_t fn, pllfloat_t dt){

     pll->v_q[0] = (pllfloat_t)(0.0);
     pll->v_q[1] = (pllfloat_t)(0.0);
  
     pll->ylf[0] = (pllfloat_t)(0.0);
     pll->ylf[1] = (pllfloat_t)(0.0);
  
     pll->fo = (pllfloat_t)(0.0);
     pll->fn = (pllfloat_t)(fn);
  
     pll->theta[0] = (pllfloat_t)(0.0);
     pll->theta[1] = (pllfloat_t)(0.0);
  
     pll->delta_t = (pllfloat_t)dt;
}
//---------------------------------------------------------------------------
void pllRun(pll_t *pll, pllfloat_t v_q){

     //
     // Update the spll_obj->v_q[0] with the grid value
     //
	pll->v_q[0] = v_q;
  
     //
     // Loop Filter
     //
	pll->ylf[0] = pll->ylf[1]
				+ (pll->lpf_coeff.b0 * pll->v_q[0])
                + (pll->lpf_coeff.b1 * pll->v_q[1]);
	pll->ylf[1] = pll->ylf[0];
	pll->v_q[1] = pll->v_q[0];
  
	pll->ylf[0] = (pll->ylf[0] > (pllfloat_t)(200.0)) ? (pllfloat_t)(200.0) : pll->ylf[0];
  
     //
     // VCO
     //
	pll->fo = pll->fn + pll->ylf[0];
  
	pll->theta[0] = pll->theta[1] + ((pll->fo * pll->delta_t) * (pllfloat_t)(2.0 * 3.1415926));

	if(pll->theta[0] > (pllfloat_t)(2.0 * 3.1415926)){
		pll->theta[0] = pll->theta[0] - (pllfloat_t)(2.0 * 3.1415926);
	}
  
     pll->theta[1] = pll->theta[0];
}
//---------------------------------------------------------------------------
void pllInitFixedPoint(pllFixedPoint_t *pll, pllint_t fn, pllint_t dt){

     pll->v_q[0] = (pllint_t)(0);
     pll->v_q[1] = (pllint_t)(0);

     pll->ylf[0] = (pllint_t)(0);
     pll->ylf[1] = (pllint_t)(0);

     pll->fo = (pllint_t)(0);
     pll->fn = (pllint_t)(fn);

     pll->theta[0] = (pllint_t)(0);
     pll->theta[1] = (pllint_t)(0);

     pll->delta_t = (pllint_t)dt;
}
//---------------------------------------------------------------------------
void pllRunFixedPoint(pllFixedPoint_t *pll, pllint_t v_q){

	//
	// Update the spll_obj->v_q[0] with the grid value
	//
	pll->v_q[0] = v_q;

	//
	// Loop Filter
	//
	pll->ylf[0] =  pll->ylf[1]
					  + fixedmul(pll->lpf_coeff.b0, pll->v_q[0])
					  + fixedmul(pll->lpf_coeff.b1, pll->v_q[1]);
	pll->ylf[1] = pll->ylf[0];
	pll->v_q[1] = pll->v_q[0];

	pll->ylf[0] = (pll->ylf[0] > fixedmathftoi(200.0f)) ? fixedmathftoi(200.0f) : pll->ylf[0];

	 //
	 // VCO
	 //
	pll->fo = pll->fn + pll->ylf[0];

	pll->theta[0] = pll->theta[1] +
			fixedmul( fixedmul(pll->fo, pll->delta_t), fixedmathftoi(2.0f * 3.1415926f));
	 if(pll->theta[0] > fixedmathftoi(2.0f * 3.1415926f))
	 {
		 pll->theta[0] = pll->theta[0] - fixedmathftoi(2.0f * 3.1415926f);
	 }

	 pll->theta[1] = pll->theta[0];
}
//---------------------------------------------------------------------------
//===========================================================================
