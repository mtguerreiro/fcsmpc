/*
 * fixedmath.h
 *
 *  Created on: 12.01.2022
 *      Author: mguerreiro
 */

#ifndef FIXEDMATH_H_
#define FIXEDMATH_H_

typedef int32_t fmint_t;
typedef int64_t fmintpr_t;

#define FIXED_MATH_Q		24
#define FIXED_MATH_FTOI		((float)(1 << FIXED_MATH_Q))
#define FIXED_MATH_ITOF		( (float) ( 1.0f / ((float)(1 << FIXED_MATH_Q)) ) )

#define fixedmul(a, b)		( (fmint_t) ( ( ((fmintpr_t)(a)) * ((fmintpr_t)(b)) ) >> FIXED_MATH_Q ) )
#define fixedadd(a, b)		( a + b )

#define fixedmathftoi(a)	( (fmint_t) (a * FIXED_MATH_FTOI) )
#define fixedmathitof(a)	( (float) (a * FIXED_MATH_ITOF) )

#endif /* FIXEDMATH_H_ */
