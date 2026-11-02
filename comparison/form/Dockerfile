# The FORM comparison image: the dev image (pqcf-dev) plus FORM 5 with FLINT 3.2.2, built on the same Debian.
#   docker build -t pqcf-dev . && docker build -t pqcf-form -f comparison/form/Dockerfile .
# FORM (GPL-3) is a separate program here: PQCF calls it as an external process (pq_form), never links it.
FROM erlang:27 AS form
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends git ca-certificates build-essential autoconf \
      automake libtool libgmp-dev libmpfr-dev zlib1g-dev >/dev/null
ARG FLINT_VERSION=v3.2.2
ARG FORM_REF=master
RUN git clone -q --depth 1 --branch "$FLINT_VERSION" https://github.com/flintlib/flint.git /flint && cd /flint \
 && ./bootstrap.sh >/dev/null 2>&1 && ./configure --prefix=/opt/form >/dev/null && make -j8 -s >/dev/null 2>&1 \
 && make install -s >/dev/null 2>&1
RUN git clone -q https://github.com/form-dev/form.git /form && cd /form && git checkout -q "$FORM_REF" \
 && autoreconf -i >/dev/null 2>&1 \
 && CPPFLAGS=-I/opt/form/include LDFLAGS="-L/opt/form/lib -Wl,-rpath,/opt/form/lib" ./configure --prefix=/opt/form > /form/configure.log 2>&1 \
 && grep -q "^ *flint" /form/configure.log && make -j8 -s >/dev/null 2>&1 && make install -s >/dev/null 2>&1 \
 && git rev-parse --short HEAD > /opt/form/FORM_COMMIT

FROM pqcf-dev
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends libgmp10 libmpfr6 zlib1g >/dev/null \
 && rm -rf /var/lib/apt/lists/*
COPY --from=form /opt/form /opt/form
ENV PATH=/opt/form/bin:$PATH
