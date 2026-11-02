FROM debian:trixie
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends git ca-certificates build-essential autoconf automake libtool libgmp-dev libmpfr-dev zlib1g-dev time python3 >/dev/null
RUN git clone -q --depth 1 --branch v3.2.2 https://github.com/flintlib/flint.git /flint && cd /flint && ./bootstrap.sh >/dev/null 2>&1 && ./configure --prefix=/usr/local >/dev/null && make -j16 -s >/dev/null 2>&1 && make install -s >/dev/null 2>&1 && ldconfig
RUN git clone -q --depth 1 https://github.com/form-dev/form.git /form && cd /form && autoreconf -i >/dev/null 2>&1 && ./configure > /configure.log 2>&1 && make -j8 -s >/dev/null 2>&1 && cp sources/form sources/tform /usr/local/bin/
