/*
 * psdtypes.h
 *
 *  Created on: 28 de dez de 2021
 *      Author: marco
 */

#ifndef PSDTYPES_H_
#define PSDTYPES_H_


typedef struct{
	float a;
	float b;
	float c;
}psdtypesABC_t;

typedef struct{
	float d;
	float q;
	float z;
}psdtypesDQ0_t;

typedef struct{
	psdtypesDQ0_t ii;
	psdtypesDQ0_t ig;
	psdtypesDQ0_t vc;
}psdtypes2LLCLDQ0Data_t;

#endif /* PSDTYPES_H_ */
