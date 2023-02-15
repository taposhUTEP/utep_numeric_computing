
/* Copyright (C) 2023 University of Texas at El Paso

   Contributed by: Christoph Lauter 
                   
                   and the 2023 class of CS4390/5390

		   Applied Numerical Computing for Multimedia
		   Applications.

   All rights reserved.

   NO LICENSE SPECIFIED.

*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "integer_ops.h"
#include "widefloat_ops.h"

/* Helper functions */

/* Tries to multiply the two size_t arguments a and b.

   If the product holds on a size_t variable, sets the 
   variable pointed to by c to that product and returns a 
   non-zero value.
   
   Otherwise, does not touch the variable pointed to by c and 
   returns zero.

   This implementation is kind of naive as it uses a division.
   If performance is an issue, try to speed it up by avoiding 
   the division while making sure that it still does the right 
   thing (which is hard to prove).

*/
static inline int __try_size_t_multiply(size_t *c, size_t a, size_t b) {
  size_t t, r, q, M;

  /* If any of the arguments a and b is zero, everthing works just fine. */
  if ((a == ((size_t) 0)) ||
      (b == ((size_t) 0))) {
    *c = a * b;
    return 1;
  }
  
  /* If both a and b are less than 2^(k/2), where k is the bitwith of 
     a size_t, a regular multiplication is enough.
  */
  M = ((size_t) 1) << (((size_t) 4) * sizeof(size_t));
  if ((a < M) && (b < M)) {
    *c = a * b;
    return 1;
  }
  
  /* Here, neither a nor b is zero. 

     We perform the multiplication, which may overflow, i.e. present
     some modulo-behavior.

  */
  t = a * b;

  /* Perform Euclidian division on t by a:

     t = a * q + r

     As we are sure that a is non-zero, we are sure
     that we will not divide by zero.

  */
  q = t / a;
  r = t % a;

  /* If the rest r is non-zero, the multiplication overflowed. */
  if (r != ((size_t) 0)) return 0;

  /* Here the rest r is zero, so we are sure that t = a * q.

     If q is different from b, the multiplication overflowed.
     Otherwise we are sure that t = a * b.

  */
  if (q != b) return 0;
  *c = t;
  return 1;
}

static inline void __m_memset(void *s, int c, size_t m, size_t n) {
  size_t p, i, min_m_n, max_m_n;
  void *curr;

  /* Easy case for 99.9999% of all cases */
  if (__try_size_t_multiply(&p, m, n)) {
    memset(s, c, p);
    return;
  }

  /* Overflow case */
  if (m < n) {
    min_m_n = m;
    max_m_n = n;
  } else {
    min_m_n = n;
    max_m_n = m;
  }
  for (i=0,curr=s;
       i<min_m_n;
       i++,curr=(void *) (((char *) curr) + max_m_n)) {
    memset(curr, c, max_m_n);
  }
}

static inline void __m_memcpy(void *dst, const void *src, size_t m, size_t n) {
  size_t p, i, min_m_n, max_m_n;
  void *curr_dst;
  const void *curr_src;

  /* Easy case for 99.9999% of all cases */
  if (__try_size_t_multiply(&p, m, n)) {
    memcpy(dst, src, p);
    return;
  }

  /* Overflow case */
  if (m < n) {
    min_m_n = m;
    max_m_n = n;
  } else {
    min_m_n = n;
    max_m_n = m;
  }
  for (i=0,curr_dst=dst,curr_src=src;
       i<min_m_n;
       i++,
	 curr_dst=(void *) (((char *) curr_dst) + max_m_n),
	 curr_src=(const void *) (((const char *) curr_src) + max_m_n)) {
    memcpy(curr_dst, curr_src, max_m_n);
  }
}

static inline void *__alloc_mem(size_t nmemb, size_t size) {
  void *ptr;

  ptr = calloc(nmemb, size);
  if (ptr == NULL) {
    fprintf(stderr, "Cannot allocate memory: %s\n", strerror(errno));
    exit(1);
  }
  
  return ptr;
}

static inline void __free_mem(void *ptr) {
  free(ptr);
}


/* Functions to (de-)allocate memory for wide floating-point
   numbers 
*/

/* Initializes the widefloat_t op to NaN, allocating memory for a
   mantissa of 64 * n - WIDEFLOAT_OVERHEAD bits of precision.
   
   Does nothing if n is zero.

*/
void widefloat_init(widefloat_t * op, size_t n) {
  if (n == ((size_t) 0)) return;
  op->fpclass = FPCLASS_NAN;
  op->sign = 0;
  op->exponent = (int32_t) 0;
  op->mantissa_size = n;
  op->mantissa = __alloc_mem(n, sizeof(*(op->mantissa)));
}


/* Deallocates the memory in the mantissa of widefloat_t op */
void widefloat_clear(widefloat_t *op) {
  op->fpclass = FPCLASS_NAN;
  op->sign = 0;
  op->exponent = (int32_t) 0;
  op->mantissa_size = (size_t) 0;
  __free_mem(op->mantissa);
  op->mantissa = NULL;
}



/* General rounding function 

   Sets the floating-point number op to the value equal to or in
   magnitude just below

   (-1)^s * 2^E * m 

   where 

   m is an integer with n digits.

   The floating-point number op must be initialized.
   
   Does nothing if op is clearly not initialized.
   
   If n is zero, sets op to zero.

   If E is too great, sets op to infinity.
   If E is too small, sets op to zero.

*/
void widefloat_set_from_scaled_integer(widefloat_t *op,
				       int s,
				       int64_t E,
				       const uint64_t *m,
				       size_t n) {
  uint64_t *t;
  int64_t EE;
  size_t q;
  uint64_t lzc;
  size_t sigma;
  int32_t expo;
  
  /* Check special cases */
  if (op->mantissa_size == ((size_t) 0)) return;
  if (op->mantissa == NULL) return;
  if (n == ((size_t) 0)) {
    op->fpclass = FP_CLASS_NUMBER;
    op->sign = !!s;
    op->exponent = (int32_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    return;
  }
  if (is_zero(m, n)) {
    op->fpclass = FP_CLASS_NUMBER;
    op->sign = !!s;
    op->exponent = (int32_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    return;
  }

  /* Here m has at least one digit and is not zero. The floating-point
     number is initialized.

     We copy m into a temporary t with size q = n + 1, so that we can
     normalize it.

  */
  q = n + ((size_t) 1);
  t = __alloc_mem(q, sizeof(*t));
  t[q - ((size_t) 1)] = (uint64_t) 0;
  __m_memcpy(t, m, n, sizeof(*t));

  /* We get a leading-zero count on t */
  lzc = leading_zeros(t, q);

  /* We shift t to the left by lzc - WIDEFLOAT_OVERHEAD */
  sigma = ((size_t) lzc) - ((size_t) WIDEFLOAT_OVERHEAD);
  shift_left(t, q, sigma);

  /* We adapt EE so that 2^EE * t = 2^E * m */
  EE = E - ((int64_t) sigma);

  /* We adapt EE to reflect a mantissa between 1 and 2 */
  EE += ((int64_t) (op->mantissa_size << 6)) -
    ((int64_t) WIDEFLOAT_OVERHEAD) -
    ((int64_t) 1);

  /* Now we check if EE holds on a 32 bit signed integer.

     If EE is greater than the greatest 32 bit signed integer, 
     we produce infinity.

     If EE is less than the smallest 32 bit signed integer, 
     we produce zero.

  */
  if (EE > ((int64_t) ((((uint64_t) 1) << 31) - ((uint64_t) 1)))) {
    /* Produce infinity */
    if (s) {
      op->fpclass = FPCLASS_NEG_INF;
    } else {
      op->fpclass = FPCLASS_POS_INF;
    }
    op->sign = !!s;
    op->exponent = (int64_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    __free_mem(t);
    return;
  }
  if (EE < ((-((int64_t) ((((uint64_t) 1) << 31) - ((uint64_t) 1)))) - ((int64_t) 1))) {
    /* Produce zero */
    op->fpclass = FP_CLASS_NUMBER;
    op->sign = !!s;
    op->exponent = (int32_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    __free_mem(t);
    return;
  }

  /* Now we know that the exponent holds on a 32bit signed integer */
  expo = (int32_t) EE;

  /* Now we can store the sign, the exponent and the mantissa in 
     the floating-point number.

     As we do round-to-zero, rounding is just a truncation.

  */
  op->fpclass = FP_CLASS_NUMBER;
  op->sign = !!s;
  op->exponent = expo;
  if (op->mantissa_size >= q) {
    /* The mantissa is longer than the temporary t */
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    __m_memcpy(&(op->mantissa[op->mantissa_size - q]), t,
	       q, sizeof(*(op->mantissa)));
  } else {
    /* The mantissa is shorter than the temporary t */
    __m_memcpy(op->mantissa, &t[q - op->mantissa_size],
	       op->mantissa_size, sizeof(*(op->mantissa)));
  }
  __free_mem(t);
}


/* Set a floating-point number from an integer

   Sets the floating-point number op to the value equal to or in
   magnitude just below

   m 

   where 

   m is an integer with n digits.

   The floating-point number op must be initialized.
   
   Does nothing if op is clearly not initialized.
   
   If n is zero, sets op to zero.

*/
void widefloat_set_from_integer(widefloat_t *op,
				const uint64_t *m,
				size_t n) {
  widefloat_set_from_scaled_integer(op, 0, (int64_t) 0, m, n);
}

