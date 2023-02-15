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
  FPCLASS_NAN      = 0,
  FPCLASS_POS_INF,
  FPCLASS_NEG_INF,
  FP_CLASS_NUMBER
} widefloatclass_t;

typedef struct {
  widefloatclass_t fpclass;
  unsigned int     sign:1;
  int32_t          exponent;
  size_t           mantissa_size;
  uint64_t         *mantissa;
} widefloat_t;

#define WIDEFLOAT_OVERHEAD  ((uint64_t) 11)

void widefloat_init(widefloat_t * op, size_t n);

void widefloat_clear(widefloat_t *op);

void widefloat_set_from_scaled_integer(widefloat_t *op,
				       int s,
				       int64_t E,
				       const uint64_t *m,
				       size_t n);

void widefloat_set_from_integer(widefloat_t *op,
				const uint64_t *m,
				size_t n);


#endif


