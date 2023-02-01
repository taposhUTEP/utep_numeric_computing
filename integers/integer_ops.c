
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

/* cout * 2^64 + s = a + b */
static inline void __halfadder(uint64_t *cout, uint64_t *s,
			       uint64_t a, uint64_t b) {
  uint64_t ccout, ss;

  ss = a + b;
  ccout = (uint64_t) (ss < a);
  *cout = ccout;
  *s = ss;
}

/* cout * 2^64 + s = a + b + cin */
static inline void __fulladder(uint64_t *cout, uint64_t *s,
			       uint64_t a, uint64_t b, uint64_t cin) {
  uint64_t cc1, cc2, ss1, ss2;

  __halfadder(&cc1, &ss1, a, b);
  __halfadder(&cc2, &ss2, ss1, cin);
  *cout = cc1 + cc2;
  *s = ss2;
}

/* s = (a + b) mod 2^(64 * k) 

   a has size m

   b has size n

   k = max(m, n)

*/
void addition(uint64_t *s,
	      const uint64_t *a, size_t m,
	      const uint64_t *b, size_t n) {
  size_t i;
  uint64_t cin, cout;

  if (m <= n) {
    /* a is shorter than b */
    cin = (uint64_t) 0;
    for (i=0;i<m;i++) {
      __fulladder(&cout, &s[i], a[i], b[i], cin);
      cin = cout;
    }
    /* Invent zeros for a */
    for (;i<n;i++) {
      __halfadder(&cout, &s[i], b[i], cin);
      cin = cout;
    }
  } else {
    /* b is shorter than a */
    cin = (uint64_t) 0;
    for (i=0;i<n;i++) {
      __fulladder(&cout, &s[i], a[i], b[i], cin);
      cin = cout;
    }
    /* Invent zeros for b */
    for (;i<m;i++) {
      __halfadder(&cout, &s[i], a[i], cin);
      cin = cout;
    }
  }
}

/* s = (a - b) mod 2^(64 * k) 

   a has size m

   b has size n

   k = max(m, n)

   Subtraction can be implemented as addition.

   a - b = a + (-b) = a + (~b + 1) = a + ~b + 1

   where ~b is a bit flip for each bit in b. 

   The extra 1 can be put into the starting carry.

*/
void subtraction(uint64_t *s,
		 const uint64_t *a, size_t m,
		 const uint64_t *b, size_t n) {
  size_t i;
  uint64_t cin, cout;

  if (m <= n) {
    /* a is shorter than b */
    cin = (uint64_t) 1;
    for (i=0;i<m;i++) {
      __fulladder(&cout, &s[i], a[i], ~b[i], cin);
      cin = cout;
    }
    /* Invent zeros for a */
    for (;i<n;i++) {
      __halfadder(&cout, &s[i], ~b[i], cin);
      cin = cout;
    }
  } else {
    /* b is shorter than a */
    cin = (uint64_t) 1;
    for (i=0;i<n;i++) {
      __fulladder(&cout, &s[i], a[i], ~b[i], cin);
      cin = cout;
    }
    /* Invent zeros for b */
    for (;i<m;i++) {
      __fulladder(&cout, &s[i], a[i], ~((uint64_t) 0), cin);
      cin = cout;
    }
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

     Finally perform a bit shift by u bits.

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

/* a gets modified such that 

   a' = floor(a / 2^k)

*/
void shift_right(uint64_t *a, size_t n, size_t k) {
  size_t w, u, i, l;
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

     Finally perform a bit shift by u bits.

     We have

     k = w * 64 + u

     This means

     u = k - 64 * w

  */
  w = k >> 6;
  u = k - (w << 6);

  /* Do the word shift, if needed */
  if (w >= ((size_t) 1)) {
    for (i=0;i<w;i++) {
      a[i] = a[i+w];
    }
    for (i=(n-w);i<n;i++) {
      a[i] = (uint64_t) 0;
    }
  }

  /* Do the bit shift, if needed */
  if (u >= ((size_t) 1)) {
    cin = (uint64_t) 0;
    for (l=n,i=(n-((uint64_t) 1));l>=1;l--,i--) {
      cout = a[i] << (64 - u);
      a[i] = (a[i] >> u) | cin;
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
  uint64_t t[n];
 
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
    addition(t, a_eight, n, a_two, n);
    addition(a, t, n, &digit, 1);
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

/* hi * 2^64 + lo = a * b */
static inline void __multiply_digits(uint64_t *hi, uint64_t *lo,
				     uint64_t a, uint64_t b) {
  uint32_t ah, al, bh, bl;
  uint64_t tah, tal, tbh, tbl;
  uint64_t hh, hl, lh, ll;
  uint32_t llh, lll, lhh, lhl, hlh, hll;
  uint64_t tllh, tlll, tlhh, tlhl, thlh, thll;
  uint32_t t1h, t1l;
  uint64_t t0, t1, t2, tt1h, tt1l, h, l;
  
  /* Cut into parts */
  ah = (uint32_t) (a >> 32);
  al = (uint32_t) a;
  bh = (uint32_t) (b >> 32);
  bl = (uint32_t) b;
  tah = (uint64_t) ah;
  tal = (uint64_t) al;
  tbh = (uint64_t) bh;
  tbl = (uint64_t) bl;
  
  /* Partial products */
  hh = tah * tbh;
  hl = tah * tbl;
  lh = tal * tbh;
  ll = tal * tbl;

  /* Cut ll, hl, lh */
  llh = (uint32_t) (ll >> 32);
  lll = (uint32_t) ll;
  hlh = (uint32_t) (hl >> 32);
  hll = (uint32_t) hl;
  lhh = (uint32_t) (lh >> 32);
  lhl = (uint32_t) lh;
  tllh = (uint64_t) llh;
  tlll = (uint64_t) lll;
  thlh = (uint64_t) hlh;
  thll = (uint64_t) hll;
  tlhh = (uint64_t) lhh;
  tlhl = (uint64_t) lhl;

  /* Add parts of same weight */
  t0 = tlll;
  t1 = tllh + thll + tlhl;
  t2 = hh + thlh + tlhh;

  /* Cut t1 into two parts */
  t1h = (uint32_t) (t1 >> 32);
  t1l = (uint32_t) t1;
  tt1h = (uint64_t) t1h;
  tt1l = (uint64_t) t1l;

  /* Add stuff together, reweighting */
  l = t0 + (tt1l << 32);
  h = t2 + tt1h;

  /* Write results back */
  *hi = h;
  *lo = l;
}

/* p = a * b

   a is on m digits
   b is on 1 digit

   p is on m+1 digits

*/
static inline void __multiplication_rectangular(uint64_t *p,
						const uint64_t *a,
						size_t m,
						uint64_t b) {
  // TODO
}

static inline unsigned int __floor_log2_size(size_t x) {
  int k;
  size_t t;
  
  /* If x is zero, return the least value we can return */
  if (x == ((size_t) 0)) return (unsigned int) 0;

  /* If x is one, return 0. */
  if (x == ((size_t) 1)) return (unsigned int) 0;

  /* Here x is at least 2

     There are all kind of smart algorithms.
     We do this with a kind of suboptimal loop.
   
  */
  for (k=((((int) sizeof(x)) * 8) - 1);
       k>=0;
       k--) {
    t = ((size_t) 1) << ((unsigned int) k);
    if ((x & t) != ((size_t) 0)) {
      return (unsigned int) k;
    }
  }
  /* Unreachable */
  return (unsigned int) 0;
}

/* r = hi * 2^(64 * s) + mi * 2^(64 * t) + lo

   where 

   hi has size m,
   mi has size n,
   lo has size p and
   r has size q.

*/
static inline  void __multiplication_helper_sum(uint64_t *r,
						size_t q,
						const uint64_t *hi,
						size_t m, size_t s,
						const uint64_t *mi,
						size_t n, size_t t,
						const uint64_t *lo,
						size_t p) {
  // TODO
}

/* p = a * b

   a is on 2^k digits
   b is on 2^k digits

   p is on 2^(k+1) digits

   with

   0 <= k <= 8 * sizeof(size_t) - 1.

*/
static inline void __multiplication_square_karatsuba(uint64_t *p,
						     const uint64_t *a,
						     const uint64_t *b,
						     unsigned int k) {
  // TODO 
}

/* Forward declaration */
static inline void __multiplication_square(uint64_t *p,
					   const uint64_t *a,
					   const uint64_t *b,
					   size_t m);

/* p = a * b

   a is on m digits
   b is on m digits

   p is on 2*m digits

   t > m

   t is a power of 2

   The function works only for t > m.

*/
static inline void __multiplication_square_aux1(uint64_t *p,
						const uint64_t *a,
						const uint64_t *b,
						size_t m,
						size_t t) {
  uint64_t aa[t];
  uint64_t bb[t];
  uint64_t r[t + t];
  
  /* Check the pre-condition */
  if (!(t > m)) return;

  /* Extend a and b to size t */
  memcpy(aa, a, m * sizeof(*aa));
  memset(&aa[m], 0, (t - m) * sizeof(*aa));
  memcpy(bb, b, m * sizeof(*bb));
  memset(&bb[m], 0, (t - m) * sizeof(*bb));

  /* Compute r = aa * bb */
  __multiplication_square(r, aa, bb, t);

  /* Compute the 2 * m last digits from r into p */
  memcpy(p, r, (m + m) * sizeof(*p));
}

/* Forward declaration */
void multiplication(uint64_t *p,
		    const uint64_t *a, size_t m,
		    const uint64_t *b, size_t n);


/* p = a * b

   a is on m digits
   b is on m digits

   p is on 2*m digits

   m > t

   t is a power of 2

   The function works only for m > t.

   On a 64 bit system, this function is essentially never needed. So
   its performance can be pretty bad.

*/
static inline void __multiplication_square_aux2(uint64_t *p,
						const uint64_t *a,
						const uint64_t *b,
						size_t m,
						size_t t) {
  uint64_t al[t];
  uint64_t ah[m - t];
  uint64_t bl[t];
  uint64_t bh[m - t];
  uint64_t hh[(m - t) + (m - t)];
  uint64_t hl[m];
  uint64_t lh[m];
  uint64_t ll[t + t];
  uint64_t hle[m + ((size_t) 1)];  /* +1 may overflow */
  uint64_t lhe[m + ((size_t) 1)];  /* +1 may overflow */
  uint64_t hllh[m + ((size_t) 1)]; /* +1 may overflow */
  
  /* Check the pre-condition */
  if (!(m > t)) return;
  
  /* Cut 

     a = ah * 2^(64 * t) + al 

     and 

     b = bh * 2^(64 * t) + bl

  */ 
  memcpy(al, a, t * sizeof(*al));
  memcpy(ah, &a[t], (m - t) * sizeof(*ah));
  memcpy(bl, b, t * sizeof(*bl));
  memcpy(bh, &b[t], (m - t) * sizeof(*bh));

  /* Perform the 4 partial products */
  multiplication(hh, ah, m - t, bh, m - t);
  multiplication(hl, ah, m - t, bl, t);
  multiplication(lh, al, t, bh, m - t);
  multiplication(ll, al, t, bl, t);

  /* Add the middle products together */
  memset(hle, 0, (m + ((size_t) 1)) * sizeof(*hle));
  memcpy(hle, hl, m * sizeof(*hle));
  memset(lhe, 0, (m + ((size_t) 1)) * sizeof(*lhe));
  memcpy(lhe, lh, m * sizeof(*lhe));
  addition(hllh, hle, m + ((size_t) 1), lhe, m + ((size_t) 1));

  /* Put everything back into the result */
  __multiplication_helper_sum(p,
			      m + m,
			      hh, (m - t) + (m - t), t + t,
			      hllh, m + ((size_t) 1), t,
			      ll, t + t);
}

/* p = a * b

   a is on m digits
   b is on m digits

   p is on 2*m digits

*/
static inline void __multiplication_square(uint64_t *p,
					   const uint64_t *a,
					   const uint64_t *b,
					   size_t m) {
  unsigned int k;
  size_t t, tt, T;

  /* If m is zero, do nothing */
  if (m == ((size_t) 0)) return;

  /* If m is one, just call the digit multiplication */
  if (m == ((size_t) 1)) {
    __multiply_digits(&p[1], &p[0], a[0], b[0]);
    return;
  }

  /* Here, m >= 2 

     Compute 

     k = floor(log2(m))

     We have 1 <= k <= 8 * sizeof(size_t) - 1.

  */
  k = __floor_log2_size(m);

  /* If m is equal to 2^k, we can call Karatsuba directly. */
  t = ((size_t) 1) << k;
  if (m == t) {
    __multiplication_square_karatsuba(p, a, b, k);
    return;
  }

  /* Here 2^k = t < m < 2^(k + 1) 

     Compute T, the largest power of 2 such that
     2 * T is representable on a size_t.

  */
  T = ((~((size_t) 0)) >> 2) + ((size_t) 1);

  /* If t is less than T, we can compute 2 * t (and 4 * t) without
     overflowing. 
  */
  if (t < T) {
    /* We can compute tt = 2 * t, which is a power of 2 and that is
       greater than m.

       We extend a and b to size tt and call ourselves.

       We do all this in a helper function.

    */
    tt = t << 1;
    __multiplication_square_aux1(p, a, b, m, tt);
    return;
  }

  /* Here m > t >= T and T is the largest power of 2 such 
     that 2 * T holds on a size_t.

     We need to cut a and b into smaller parts. We do so in 
     a helper function.

  */
  __multiplication_square_aux2(p, a, b, m, T);
}

/* p = a * b 

   Works only for:

   +  a on m digits
   +  b on n digits 
   +  2 <= m
   +  m < n

*/
static inline void __multiplication_aux(uint64_t *p,
					const uint64_t *a,
					size_t m,
					const uint64_t *b,
					size_t n) {
  uint64_t t[n];
  uint64_t r[n+n]; /* n+n may overflow */
  
  /* Handle preconditions */
  if (!(((size_t) 2) <= m)) return;
  if (!(m < n)) return;

  /* Copy a into t */
  memcpy(t, a, m * sizeof(*t));

  /* Set upper part of t to zero */
  memset(&t[m], 0, (n - m) * sizeof(*t));

  /* Now t and b have the same size n */
  __multiplication_square(r, t, b, n);

  /* Copy low part of r into p */
  memcpy(p, r, (m + n) * sizeof(*p));
}

/* p = a * b

   a is on m "digits"
   b is on n "digits"

   p must have m + n "digits"

*/
void multiplication(uint64_t *p,
		    const uint64_t *a, size_t m,
		    const uint64_t *b, size_t n) {
  
  /* If one of the sizes is zero, we do nothing */
  if (m == ((size_t) 0)) return;
  if (n == ((size_t) 0)) return;
  
  /* If m > n, we flip the arguments */
  if (m > n) {
    multiplication(p, b, n, a, m);
    return;
  }

  /* Here, m <= n

     If m is 1 or n is 1, call a rectangular 
     schoolbook multiplication.

  */
  if (m == ((size_t) 1)) {
    __multiplication_rectangular(p, b, n, a[0]);
    return;
  }
  if (n == ((size_t) 1)) {
    __multiplication_rectangular(p, a, m, b[0]);
    return;
  }

  /* Here, 2 <= m <= n.

     If m = n, we call a square multiplication.

  */
  if (m == n) {
    __multiplication_square(p, a, b, m);
    return;
  }

  /* Here, 2 <= m < n.

     We extend a to the size of n.

     We do all that in a helper function.

  */
  __multiplication_aux(p, a, m, b, n);
}

