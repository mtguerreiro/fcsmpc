/*
* tlvsi_config.h
*
*  Created on: 23 de ago de 2022
*      Author: marco
*/

#ifndef TLVSI_CONFIG_H_
#define TLVSI_CONFIG_H_

//===========================================================================
/*------------------------------- Includes --------------------------------*/
//===========================================================================
#include <stdint.h>

//===========================================================================

//===========================================================================
/*------------------------------ Definitions ------------------------------*/
//===========================================================================
#define TLVSI_CONFIG_SCALE		200.0f

#define TLVSI_CONFIG_Li			((float)0.0034)
#define TLVSI_CONFIG_Lg			((float)0.0018)
#define TLVSI_CONFIG_Cf			((float)2e-05)

#define TLVSI_CONFIG_V_dc		((float)(650.0f / TLVSI_CONFIG_SCALE))

#define TLVSI_CONFIG_w			((float)314.1592653589793)

#define TLVSI_CONFIG_ts			((float)1.0f/40000.0f)

#define TLVSI_CONFIG_w_ii		((float)(1.0f / TLVSI_CONFIG_SCALE))
#define TLVSI_CONFIG_w_ig		((float)(625.0f / TLVSI_CONFIG_SCALE))
#define TLVSI_CONFIG_w_vc		((float)(0.48999999999999994f / TLVSI_CONFIG_SCALE))

//===========================================================================

#endif /* TLVSI_CONFIG_H_ */