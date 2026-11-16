// one integral per input line; prints Re Im of the finite part (mu^2 = 1) with 17 significant digits
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "qcdloop/qcdloop.h"
using std::vector;
int main() {
  ql::Triangle<std::complex<double>, double, double> tr;
  ql::Box<std::complex<double>, double, double> bx;
  vector<std::complex<double>> res(3);
  std::string kind;
  while (std::cin >> kind) {
    if (kind == "C") {
      vector<double> p(3), m(3);
      for (auto &x : p) std::cin >> x;
      for (auto &x : m) std::cin >> x;
      tr.integral(res, 1.0, m, p);
    } else {
      vector<double> p(6), m(4);
      for (auto &x : p) std::cin >> x;
      for (auto &x : m) std::cin >> x;
      bx.integral(res, 1.0, m, p);
    }
    std::printf("%.17g %.17g\n", res[0].real(), res[0].imag());
  }
  return 0;
}
