# QCDLoop 2.1.0 (Carrazza, Ellis and Zanderighi, Comput. Phys. Commun. 209, 134 (2016), arXiv:1605.03181; GPL-3) for
# the comparison of PQCF's three- and four-point functions (plan "theoretical physics instruments", F3). The release
# tarball is checked by its SHA-256. ql_eval reads lines "C p1 p2 p3 m1 m2 m3" / "D p1 p2 p3 p4 s12 s23 m1 m2 m3 m4"
# (squared, LoopTools' order) and prints the finite part at mu^2 = 1 to 17 digits.
#   docker build -t pqcf-qcdloop comparison/oneloop/qcdloop
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends g++ cmake make curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
ARG QL_SHA=75ac2d00231714210859d39c4beebce14f183707afe47f4550f33ef9f348b3c1
RUN curl -sL -o /tmp/ql.tgz https://github.com/scarrazza/qcdloop/archive/refs/tags/2.1.0.tar.gz \
 && echo "$QL_SHA  /tmp/ql.tgz" | sha256sum -c - && tar xzf /tmp/ql.tgz -C /opt && cd /opt/qcdloop-2.1.0 \
 && mkdir build && cd build && cmake -DENABLE_FORTRAN_WRAPPER=OFF -DCMAKE_INSTALL_PREFIX=/opt/ql .. >/dev/null \
 && make -j8 >/dev/null && make install >/dev/null
COPY ql_eval.cc /opt/ql_eval.cc
RUN g++ -std=c++11 -O2 /opt/ql_eval.cc -I/opt/ql/include -L/opt/ql/lib -lqcdloop -o /usr/local/bin/ql_eval
ENV LD_LIBRARY_PATH=/opt/ql/lib
