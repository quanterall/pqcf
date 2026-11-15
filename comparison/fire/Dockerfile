# FIRE 6.5.2 with the FLINT simplifier (no Fermat), for the comparison in decision record R1. On aarch64 it needs
# -fsigned-char (FIRE's t_index is a plain char); the PARI submodule is skipped (FLINT doesn't need it).
# FIRE is pinned at the commit our published runs used (its submodules follow FIRE's own pins).
#   docker build -t pqcf-fire2 comparison/fire          (the name bench/ibp.sh expects)
#   docker run --rm -v $PWD/comparison/fire:/cfg pqcf-fire2 sh -c 'cp /cfg/doublebox1.config examples/ && bin/FIRE6 -c examples/doublebox1'
FROM debian:bookworm
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    g++ make cmake autoconf automake libtool pkg-config git ca-certificates python3 wget curl patch xz-utils m4 bzip2 \
    libflint-dev libgmp-dev libmpfr-dev zlib1g-dev libzstd-dev && rm -rf /var/lib/apt/lists/*
ARG FIRE_REF=51bd31ac3496bc0e88399ee32e4429a500fa796d
RUN git clone https://gitlab.com/feynmanintegrals/fire.git /opt/fire && cd /opt/fire && git checkout -q "$FIRE_REF"
WORKDIR /opt/fire/FIRE6
RUN cd .. && git submodule update --init FIRE6/extra/fuel && cd FIRE6/extra/fuel && git submodule update --init --depth 1 extra/flint
RUN ./configure --enable_flint && \
    grep -rl "git submodule update" --include="Makefile*" . | xargs sed -i 's/git submodule update --init --recursive/true/' && \
    grep -rl -- "-O3 -Wall" --include="Makefile*" . | xargs sed -i 's/-O3 -Wall/-O3 -fsigned-char -Wall/'
RUN make dep -j16 && make -j16
