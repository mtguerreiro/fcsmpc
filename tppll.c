/*
 * tppll.c
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco.
 *
 */

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include "tppll.h"

#include <math.h>
//===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
//---------------------------------------------------------------------------
void tppllInit(float32_t grid_freq, float32_t delta_t, SPLL_3PH_SRF *spll_obj){

     spll_obj->v_q[0] = (float32_t)(0.0);
     spll_obj->v_q[1] = (float32_t)(0.0);
  
     spll_obj->ylf[0] = (float32_t)(0.0);
     spll_obj->ylf[1] = (float32_t)(0.0);
  
     spll_obj->fo = (float32_t)(0.0);
     spll_obj->fn = (float32_t)(grid_freq);
  
     spll_obj->theta[0] = (float32_t)(0.0);
     spll_obj->theta[1] = (float32_t)(0.0);
  
     spll_obj->delta_t = (float32_t)delta_t;
}
//---------------------------------------------------------------------------
//void tppllRun(float32_t v_q, SPLL_3PH_SRF *spll_obj){
//
//     //
//     // Update the spll_obj->v_q[0] with the grid value
//     //
//     spll_obj->v_q[0] = v_q;
//
//     //
//     // Loop Filter
//     //
//     spll_obj->ylf[0] =  spll_obj->ylf[1]
//                      + (spll_obj->lpf_coeff.b0 * spll_obj->v_q[0])
//                      + (spll_obj->lpf_coeff.b1 * spll_obj->v_q[1]);
//     spll_obj->ylf[1] = spll_obj->ylf[0];
//     spll_obj->v_q[1] = spll_obj->v_q[0];
//
//     spll_obj->ylf[0] = (spll_obj->ylf[0] > (float32_t)(200.0)) ?
//                                 (float32_t)(200.0) : spll_obj->ylf[0];
//
//     //
//     // VCO
//     //
//     spll_obj->fo = spll_obj->fn + spll_obj->ylf[0];
//
//     spll_obj->theta[0] = spll_obj->theta[1] +
//                          ((spll_obj->fo * spll_obj->delta_t) *
//                           (float32_t)(2.0 * 3.1415926));
//     if(spll_obj->theta[0] > (float32_t)(2.0 * 3.1415926))
//     {
//         spll_obj->theta[0] = spll_obj->theta[0] - (float32_t)(2.0 * 3.1415926);
//     }
//
//     spll_obj->theta[1] = spll_obj->theta[0];
//}
//---------------------------------------------------------------------------
void tppllRunInt(fmint_t v_q, SPLL_3PH_SRF_INT *spll_obj){

     //
     // Update the spll_obj->v_q[0] with the grid value
     //
     spll_obj->v_q[0] = v_q;

     //
     // Loop Filter
     //
     spll_obj->ylf[0] =  spll_obj->ylf[1]
                      + fixedmul(spll_obj->lpf_coeff.b0, spll_obj->v_q[0])
                      + fixedmul(spll_obj->lpf_coeff.b1, spll_obj->v_q[1]);
     spll_obj->ylf[1] = spll_obj->ylf[0];
     spll_obj->v_q[1] = spll_obj->v_q[0];

     spll_obj->ylf[0] = (spll_obj->ylf[0] > fixedmathftoi(200.0f)) ?
    		 	 	 	 	 	 fixedmathftoi(200.0f) : spll_obj->ylf[0];

     //
     // VCO
     //
     spll_obj->fo = spll_obj->fn + spll_obj->ylf[0];

     spll_obj->theta[0] = spll_obj->theta[1] +
                          fixedmul( fixedmul(spll_obj->fo, spll_obj->delta_t),
                        		fixedmathftoi(2.0f * 3.1415926f));
     if(spll_obj->theta[0] > fixedmathftoi(2.0f * 3.1415926f))
     {
         spll_obj->theta[0] = spll_obj->theta[0] - fixedmathftoi(2.0f * 3.1415926f);
     }

     spll_obj->theta[1] = spll_obj->theta[0];
}
//---------------------------------------------------------------------------
//===========================================================================
