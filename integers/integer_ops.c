
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

/* cout * 2^64 + s = a + b + cin */
static inline void __fulladder(uint64_t *cout, uint64_t *s,
			       uint64_t a, uint64_t b, uint64_t cin) {
  uint64_t ccout, ss;

  ccout = (uint64_t) 0;
  ss = a + b;
  ccout += (uint64_t) (ss < a);
  ss += cin;
  ccout += (uint64_t) (ss < cin);
  *cout = ccout;
  *s = ss;
}

/* s = (a + b) mod 2^(64 * n) */
void addition(uint64_t *s, const uint64_t *a, const uint64_t *b, size_t n) {
  size_t i;
  uint64_t cin, cout;

  cin = (uint64_t) 0;
  for (i=0;i<n;i++) {
    __fulladder(&cout, &s[i], a[i], b[i], cin);
    cin = cout;
  }
}

/* Modifies a of size n s.t.

   a' = (a + b) mod 2^(64 * n)

*/
static inline void __add_uint64(uint64_t *a, size_t n, uint64_t b) {
  size_t i;
  uint64_t cin, cout;

  if (n == ((size_t) 0)) return;
  cin = (uint64_t) 0;
  __fulladder(&cout, &a[0], a[0], b, cin);
  for (i=1;i<n;i++) {
    __fulladder(&cout, &a[i], a[i], 0, cin);
    cin = cout;
  }
}

/* a gets modified such that 

   a' = (a * 2^k) mod 2^(64 * n)

*/
void shift_left(uint64_t *a, size_t n, size_t k) {
  size_t w, u, i;
  uint64_t cin, cout;
  
  /* Shift by 0 bits => nothing to do */
  if (k == ((size_t) 0)) return;

  /* Shift by at least n * 64 bits => set 
     result to zero.
  */
  if (k >= (n * ((size_t) 64))) {
    /* Helper function to set memory to zero */
    memset(a, 0, n * sizeof(*a));
    return;
  }

  /* Here we have 1 <= k <= 64 * n - 1 

     First "cut" k into w and u such that

     k = w * 64 + u

     Then perform a word shift by w words of
     64 bits.

     Finally performa a bit shift by u bits.

     We have

     k = w * 64 + u

     This means

     u = k - 64 * w

  */
  w = k >> 6;
  u = k - (w << 6);

  /* Do the word shift, if needed */
  if (w >= ((size_t) 1)) {
    for (i=(n-((size_t) 1));i>=w;i--) {
      a[i] = a[i-w];
    }
    for (i=0;i<w;i++) {
      a[i] = (uint64_t) 0;
    }
  }

  /* Do the bit shift, if needed */
  if (u >= ((size_t) 1)) {
    cin = (uint64_t) 0;
    for (i=0;i<n;i++) {
      cout = a[i] >> (64 - u);
      a[i] = (a[i] << u) | cin;
      cin = cout;
    }
  }
}

/* a becomes what is in the decimal string mod 2^(64 * n)

   Returns 0 if success
   Returns -1 if failure  (someone tries to convert "cheese")

*/
int convert_from_decimal_string(uint64_t *a, size_t n,
				const char *str) {
  const char *curr;
  uint64_t digit;
  uint64_t a_eight[n];
  uint64_t a_two[n];
 
  /* Do nothing for the empty string */
  if (str[0] == '\0') return -1;
  
  /* Set a to zero */
  memset(a, 0, n * sizeof(*a));

  /* Loop over the string */
  for (curr=str; *curr!='\0'; curr++) {
    if (!(('0' <= *curr) &&
	  (*curr <= '9'))) return -1;
    digit = (uint64_t) (((int) *curr) - ((int) '0'));

    /* Multiply a by 10 and add in the digit 

       a * 10 = a * 8 + a * 2

    */
    memcpy(a_eight, a, n * sizeof(*a));
    memcpy(a_two, a, n * sizeof(*a));
    shift_left(a_eight, n, 3);
    shift_left(a_two, n, 1);
    addition(a, a_eight, a_two, n);
    __add_uint64(a, n, digit);
  }
  return 0;
}

/* Returns  -1  if a < b
             0  if a = b
             1  if a > b
*/
int comparison(const uint64_t *a, const uint64_t *b, size_t n) {
  size_t i, k;

  for (i=(n-((size_t) 1)),k=n;k>((size_t) 0);i--,k--) {
    if (a[i] < b[i]) return -1;
    if (a[i] > b[i]) return 1;
  }
  return 0;
}


