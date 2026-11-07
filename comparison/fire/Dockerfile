# FIRE 6.5.2 with the FLINT simplifier (no Fermat), for the comparison in decision record R1. On aarch64 it needs
# -fsigned-char (FIRE's t_index is a plain char); the PARI submodule is skipped (FLINT doesn't need it).
#   docker build -t pqcf-fire comparison/fire
#   docker run --rm -v $PWD/comparison/fire:/cfg pqcf-fire sh -c 'cp /cfg/doublebox1.config examples/ && bin/FIRE6 -c examples/doublebox1'
FROM debian:bookworm
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    g++ make cmake autoconf automake libtool pkg-config git ca-certificates python3 wget curl patch xz-utils m4 bzip2 \
    libflint-dev libgmp-dev libmpfr-dev zlib1g-dev libzstd-dev && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 https://gitlab.com/feynmanintegrals/fire.git /opt/fire
WORKDIR /opt/fire/FIRE6
RUN cd .. && git submodule update --init FIRE6/extra/fuel && cd FIRE6/extra/fuel && git submodule update --init --depth 1 extra/flint
RUN ./configure --enable_flint && \
    grep -rl "git submodule update" --include="Makefile*" . | xargs sed -i 's/git submodule update --init --recursive/true/' && \
    grep -rl -- "-O3 -Wall" --include="Makefile*" . | xargs sed -i 's/-O3 -Wall/-O3 -fsigned-char -Wall/'
RUN make dep -j16 && make -j16
