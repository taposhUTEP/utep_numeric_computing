
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

