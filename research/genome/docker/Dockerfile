# The genome-research baselines: NumPy, SciPy, Numba (simulated annealing, parallel tempering), NetworkX, and
# OR-Tools' CP-SAT as the exact solver, and scikit-learn (G3).
#   docker build -t pqcf-genome research/genome/docker
FROM python:3.11-slim
RUN pip install --no-cache-dir numpy==2.3.4 scipy==1.16.3 numba==0.62.1 networkx==3.5 ortools==9.14.6206 scikit-learn==1.7.2
ENV NUMBA_CACHE_DIR=/tmp/numba
