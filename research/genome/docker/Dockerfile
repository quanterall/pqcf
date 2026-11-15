# The genome-research baselines: NumPy, SciPy, Numba (simulated annealing, parallel tempering), NetworkX, and
# OR-Tools' CP-SAT as the exact solver, and scikit-learn (G3).
#   docker build -t pqcf-genome research/genome/docker
FROM python:3.11-slim@sha256:e88e9763f943ec1834f992a4b51e0f24500486803e8bc534e5767af9ea65f6ce
RUN pip install --no-cache-dir numpy==2.3.4 scipy==1.16.3 numba==0.62.1 networkx==3.5 ortools==9.14.6206 scikit-learn==1.7.2
ENV NUMBA_CACHE_DIR=/tmp/numba
