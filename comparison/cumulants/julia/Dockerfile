# QuantumCumulants.jl (Plankensteiner, Hotter, Ritsch, Quantum 6, 617 (2022), arXiv:2105.01657) for the comparison of
# PQCF's cumulant expansions (quantum optics plan, Q4). Versions pinned by Project.toml/Manifest written at build time.
#   docker build -t pqcf-cumulants comparison/cumulants/julia
FROM julia:1.11
ENV JULIA_PROJECT=/opt/qc JULIA_DEPOT_PATH=/opt/julia
RUN mkdir -p /opt/qc && julia -e 'using Pkg; Pkg.add([PackageSpec(name="QuantumCumulants", version="0.3"), \
        PackageSpec(name="ModelingToolkit"), PackageSpec(name="OrdinaryDiffEq"), PackageSpec(name="Symbolics")]); \
        Pkg.precompile(); using QuantumCumulants; println(pkgversion(QuantumCumulants))'
RUN julia -e 'using Pkg; Pkg.status()'
