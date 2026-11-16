# CRunDec 3.1.2 (Herren and Steinhauser, Comput. Phys. Commun. 224, 333 (2018), arXiv:1703.03751; MIT licence) at its
# authors' repository, commit 31e0c7f1085c0bbf0aa4340c0cd8dd29aae4af65, the tarball checked by SHA-256: the reference
# for PQCF's running and decoupling of alpha_s and quark masses (plan "theoretical physics instruments", QCD1).
# crd_eval reads lines and prints one number each to 17 digits:
#   A as0 mu0 mu1 nf loops        alpha_s(mu1)                 AlphasExact
#   M m0 as0 as1 nf loops         m(mu1) from m(mu0)           mMS2mMS (RunDec's expanded formula)
#   D as mh mu nl loops           alpha_s^(nl) from ^(nl+1)    DecAsDownMS
#   U as mh mu nl loops           alpha_s^(nl+1) from ^(nl)    DecAsUpMS
#   X m0 as0 mu0 mu1 nf loops     m(mu1) alpha_s(mu1), exactly  AsmMSrunexact (two numbers)
#   docker build -t pqcf-crundec comparison/qcd/crundec
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends g++ curl ca-certificates && rm -rf /var/lib/apt/lists/*
ARG CRD=31e0c7f1085c0bbf0aa4340c0cd8dd29aae4af65
ARG CRD_SHA=2c20263b5a37af1ae1024be67159efa941ad4f330766db33de0301d654688447
RUN curl -sL -o /tmp/crd.tgz https://github.com/Herren/CRunDec/archive/$CRD.tar.gz \
 && echo "$CRD_SHA  /tmp/crd.tgz" | sha256sum -c - && tar xzf /tmp/crd.tgz -C /opt && mv /opt/CRunDec-$CRD /opt/crd
COPY crd_eval.cpp /opt/crd/crd_eval.cpp
RUN cd /opt/crd && g++ -std=c++11 -O2 -I. -o /usr/local/bin/crd_eval crd_eval.cpp CRunDec.cpp
