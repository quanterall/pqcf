# The oracles of FS 12.1, pinned; used offline to generate test/fixtures, never at run time.
#   docker build -t pqcf-oracle test/oracle && docker run --rm -v "$PWD":/src -w /src pqcf-oracle python test/oracle/<script>.py
FROM python:3.11-slim
# Piquasso has no aarch64 wheel: its extension is built here
RUN apt-get update && apt-get install -y --no-install-recommends g++ make && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt
ENV NUMBA_CACHE_DIR=/tmp/numba
