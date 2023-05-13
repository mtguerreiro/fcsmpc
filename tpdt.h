/*
 * tpdt.h
 *
 *	Three-phase data types.
 *
 *  Created on: 28 de dez de 2021
 *      Author: marco
 */

#ifndef TPDT_H_
#define TPDT_H_

typedef float tpdtfloat_t;
typedef int32_t tpdtint_t;

typedef struct{
	tpdtfloat_t a;
	tpdtfloat_t b;
	tpdtfloat_t c;
}tpdtAbc_t;

typedef struct{
	tpdtint_t a;
	tpdtint_t b;
	tpdtint_t c;
}tpdtAbcFixedPoint_t;

typedef struct{
	tpdtfloat_t d;
	tpdtfloat_t q;
	tpdtfloat_t z;
}tpdtDq0_t;

typedef struct{
	tpdtint_t d;
	tpdtint_t q;
	tpdtint_t z;
}tpdtDq0FixedPoint_t;

typedef struct{
	tpdtfloat_t ab;
	tpdtfloat_t ac;
	tpdtfloat_t bc;
}tpdtAbcLine_t;

typedef struct{
	tpdtint_t ab;
	tpdtint_t ac;
	tpdtint_t bc;
}tpdtAbcLineFixedPoint_t;

#endif /* TPDT_H_ */
