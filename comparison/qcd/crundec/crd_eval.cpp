// one computation per input line (see the Dockerfile); prints the result with 17 significant digits
#include <cstdio>
#include <iostream>
#include <string>
#include "CRunDec.h"
int main() {
  CRunDec crd;
  std::string k;
  double a, b, c;
  int n, l;
  while (std::cin >> k) {
    double r;
    if (k == "X") {
      // m and alpha_s at mu1 together, exactly: X m0 as0 mu0 mu1 nf loops; prints m(mu1) then alpha_s(mu1)
      double m0, as0, mu0, mu1;
      std::cin >> m0 >> as0 >> mu0 >> mu1 >> n >> l;
      AsmMS x = crd.AsmMSrunexact(m0, as0, mu0, mu1, n, l);
      std::printf("%.17g %.17g\n", x.mMSexact, x.Asexact);
      continue;
    }
    std::cin >> a >> b >> c >> n >> l;
    if (k == "A") r = crd.AlphasExact(a, b, c, n, l);
    else if (k == "M") r = crd.mMS2mMS(a, b, c, n, l);
    else if (k == "D") r = crd.DecAsDownMS(a, b, c, n, l);
    else r = crd.DecAsUpMS(a, b, c, n, l);
    std::printf("%.17g\n", r);
  }
  return 0;
}
