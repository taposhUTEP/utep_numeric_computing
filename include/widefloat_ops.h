/* Copyright (C) 2023 University of Texas at El Paso

   Contributed by: Christoph Lauter 
                   
                   and the 2023 class of CS4390/5390

		   Applied Numerical Computing for Multimedia
		   Applications.

   All rights reserved.

   NO LICENSE SPECIFIED.

*/

#ifndef WIDEFLOAT_OPS_H
#define WIDEFLOAT_OPS_H

#include <stdint.h>

typedef enum {
  FPCLASS_NAN,
  FPCLASS_POS_INF,
  FPCLASS_NEG_INF,
  FP_CLASS_NUMBER
} widefpclass_t;

typedef struct {
  widefpclass_t fpclass;
  int32_t       exponent;
  size_t        mantissa_size;
  uint64_t      *mantissa;
} widefp_t;




#endif


