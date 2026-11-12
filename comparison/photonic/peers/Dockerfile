# The photonic peers, pinned: The Walrus, Piquasso (its C++ kernels), Perceval, DeepQuantum (PyTorch, CPU).
#   docker build -t pqcf-photonic-peers comparison/photonic/peers
FROM python:3.11-slim
# Piquasso has no aarch64 wheel: its extension is built here
RUN apt-get update && apt-get install -y --no-install-recommends g++ make cmake && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu
# Piquasso 8.0.2 keeps a sign in a plain char (src/permanent.cpp: char parity = ... ? 1 : -1), and char is unsigned on
# aarch64: its permanent is wrong there (2 x 2 [[1, 2], [3, 4]] gives 522, not 10). Built with signed chars, as on x86.
RUN CFLAGS=-fsigned-char CXXFLAGS=-fsigned-char pip install --no-cache-dir --no-binary piquasso \
        thewalrus==0.22.0 piquasso==8.0.2 perceval-quandela==1.3.1 deepquantum==4.5.0
ENV NUMBA_CACHE_DIR=/tmp/numba
