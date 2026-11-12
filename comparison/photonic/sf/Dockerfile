# Strawberry Fields 0.23.0, the last release (the repository is archived); it needs NumPy 1.x.
#   docker build -t pqcf-photonic-sf comparison/photonic/sf
FROM python:3.10-slim
RUN pip install --no-cache-dir "numpy==1.26.4" "scipy==1.11.4" "thewalrus==0.21.0" "networkx<3.3" strawberryfields==0.23.0
ENV NUMBA_CACHE_DIR=/tmp/numba
