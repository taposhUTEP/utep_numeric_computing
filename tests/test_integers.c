
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "utepnum.h"


static int convert_to_size(size_t *res, const char *str) {
  unsigned long long int t, ttt;
  size_t tt;
  char *end;

  if (str == NULL) return -1;
  if (*str == '\0') return -1;
  t = strtoull(str, &end, 0);
  if (*end != '\0') return -1;
  tt = (size_t) t;
  ttt = (unsigned long long int) tt;
  if (ttt != t) return -1;
  *res = tt;
  return 0;
}

static void print_array(const char *str, const uint64_t *a, size_t n) {
  char *sep = "";
  size_t i, k;
  
  printf("%s[", str);
  if (n == ((size_t) 0)) {
    printf("]\n");
    return;
  }

  for (k=n,i=n-((size_t) 1);k>0;k--,i--) {
    printf("%s%llu", sep, (unsigned long long int) a[i]);
    sep = ",";
  }
  printf("]\n");
}

int test_integers(size_t m, size_t n, const char *str1, const char *str2) {
  size_t q = (m > n) ? m : n;
  uint64_t a[m];
  uint64_t b[n];
  uint64_t c[q];

  /* Convert the two strings str1 and str2 */
  if (convert_from_decimal_string(a, m, str1) < 0) return -1;
  if (convert_from_decimal_string(b, n, str2) < 0) return -1;
  
  printf("Testing Addition:\n");
  /* Display the two arrays */
  print_array("a = ", a, m);
  print_array("b = ", b, n);

  /* Call addition */
  addition(c, a, m, b, n);

  /* Display the addition result */
  print_array("c = ", c, q);
  
  /* TODO */
  printf("Testing Multiplication:\n");
  /* Display the two arrays */
  print_array("a = ", a, m);
  print_array("b = ", b, n);

  /* Call multiplication */
  multiplication(c, a, m, b, n);

  /* Display the multiplication result */
  print_array("c = ", c, q);

  printf("Testing Subtraction:\n");
  /* Display the two arrays */
  print_array("a = ", a, m);
  print_array("b = ", b, n);

  /* Call Subtraction */
  subtraction(c, a, m, b, n);

  /* Display the Subtraction result */
  print_array("c = ", c, q);

  printf("Testing comparision:\n");
  /* Display the two arrays */
  print_array("a = ", a, m);
  print_array("b = ", b, n);

  /* Call comparision and display the comparision result*/
  printf("comparison result: %d\n", comparison(a, b, q));

  /* Signal success */
  return 0;
}

int main(int argc, char **argv) {
  size_t m, n;
  
  /* Check if we have at least 5 arguments */
  if (argc < 5) return 1;

  /* Convert the first two arguments to size_t */
  if (convert_to_size(&m, argv[1]) < 0) return 1;
  if (convert_to_size(&n, argv[2]) < 0) return 1;

  /* Check that none of m or n is zero */
  if (m == ((size_t) 0)) return 1;
  if (n == ((size_t) 0)) return 1;
  
  /* Run the actual test function */
  if (test_integers(m, n, argv[3], argv[4]) < 0) return 1;

  /* Signal success */
  return 0;
}
