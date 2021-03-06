/*
 * tptransforms.h
 *
 *  Created on: 25 de dez de 2021
 *      Author: marco
 *
 * The source code for the transform was extracted (and slightly modified)
 * from TI's DigitalPower library.
 */

#ifndef TPTRANSFORMS_H_
#define TPTRANSFORMS_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

//===========================================================================

////===========================================================================
///*------------------------------ Definitions ------------------------------*/
////===========================================================================
//typedef struct{
//	float a;
//	float b;
//	float c;
//
//	float d;
//	float q;
//	float z;
//}tptransformsABCDQ0_t;
////===========================================================================

//===========================================================================
/*------------------------------- Functions -------------------------------*/
//===========================================================================
void tptransformsABCDQ0(float *abc, float *dq0, float theta);
//===========================================================================

#endif /* TPTRANSFORMS_H_ */
