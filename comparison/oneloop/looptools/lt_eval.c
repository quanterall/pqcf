#include <stdio.h>
#include <string.h>
#include "clooptools.h"
/* one integral per input line; prints Re Im with 17 significant digits (for B: B0 B1 B00 B11 DB0) */
int main(void) {
  char kind[4];
  double v[10];
  ltini();
  while (scanf("%3s", kind) == 1) {
    ComplexType r;
    if (!strcmp(kind, "A")) { scanf("%lf", &v[0]); r = A0(v[0]); }
    else if (!strcmp(kind, "B")) {
      /* B0 B1 B00 B11 DB0 on one line */
      scanf("%lf %lf %lf", &v[0], &v[1], &v[2]);
      ComplexType b[5] = {B0(v[0], v[1], v[2]), B1(v[0], v[1], v[2]), B00(v[0], v[1], v[2]), B11(v[0], v[1], v[2]),
                          DB0(v[0], v[1], v[2])};
      for (int i = 0; i < 5; i++) printf("%.17g %.17g%s", Re(b[i]), Im(b[i]), i < 4 ? " " : "\n");
      continue;
    }
    else if (!strcmp(kind, "C")) {
      /* C0nocache: res[0] the finite part (res[1], res[2] the 1/eps and 1/eps^2 coefficients) */
      ComplexType res[3];
      for (int i = 0; i < 6; i++) scanf("%lf", &v[i]);
      C0nocache(res, v[0], v[1], v[2], v[3], v[4], v[5]);
      r = res[0];
    }
    else { for (int i = 0; i < 10; i++) scanf("%lf", &v[i]); r = D0(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9]); }
    printf("%.17g %.17g\n", Re(r), Im(r));
  }
  ltexi();
  return 0;
}
