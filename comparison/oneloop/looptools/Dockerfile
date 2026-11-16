# LoopTools 2.16 (Hahn and Perez-Victoria, Comput. Phys. Commun. 118, 153 (1999), arXiv:hep-ph/9807565; LGPL-3) for
# the comparison of PQCF's one-loop integrals (plan "theoretical physics instruments", F3). The source tarball is
# checked by its SHA-256. lt_eval reads lines "A m2 | B p2 m1 m2 | C p1 p2 p3 m1 m2 m3 | D ..." and prints the finite
# part (Delta = 0, mu^2 = 1, LoopTools' defaults) to 17 digits.
#   docker build -t pqcf-looptools comparison/oneloop/looptools
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends gfortran gcc g++ make curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
ARG LT_SHA=4113467575db3a14405d62d9e516b3b90410b73ea8d20c8eb8d70a30fc5cc9cb
RUN curl -sL -o /tmp/lt.tgz https://feynarts.de/looptools/LoopTools-2.16.tar.gz \
 && echo "$LT_SHA  /tmp/lt.tgz" | sha256sum -c - && tar xzf /tmp/lt.tgz -C /opt && cd /opt/LoopTools-2.16 \
 && ./configure --prefix=/opt/lt >/dev/null && make -s >/dev/null 2>&1 && make -s install >/dev/null 2>&1
COPY lt_eval.c /opt/lt_eval.c
RUN cd /opt && g++ -x c++ lt_eval.c -I/opt/lt/include -L/opt/lt/lib64 -L/opt/lt/lib -looptools -lgfortran -o /usr/local/bin/lt_eval
