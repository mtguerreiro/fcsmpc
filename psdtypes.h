/*
 * psdtypes.h
 *
 *  Created on: 28 de dez de 2021
 *      Author: marco
 */

#ifndef PSDTYPES_H_
#define PSDTYPES_H_

#include "fixedmath.h"

typedef struct{
	float a;
	float b;
	float c;
}psdtypesABC_t;

typedef struct{
	fmint_t a;
	fmint_t b;
	fmint_t c;
}psdtypesABCint_t;

typedef struct{
	float d;
	float q;
	float z;
}psdtypesDQ0_t;

typedef struct{
	fmint_t d;
	fmint_t q;
	fmint_t z;
}psdtypesDQ0int_t;

//typedef struct{
//	int32_t a;
//	int32_t b;
//	int32_t c;
//}psdtypesABCint_t;
//
//typedef struct{
//	float d;
//	float q;
//	float z;
//}psdtypesDQ0_t;
//
//typedef struct{
//	int32_t d;
//	int32_t q;
//	int32_t z;
//}psdtypesDQ0int_t;

typedef struct{
	psdtypesDQ0_t ii;
	psdtypesDQ0_t ig;
	psdtypesDQ0_t vc;
}psdtypes2LLCLDQ0Data_t;

#endif /* PSDTYPES_H_ */
