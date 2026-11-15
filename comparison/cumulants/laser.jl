# QuantumCumulants.jl's single-atom laser (its documentation's first example; arXiv:2105.01657), order 2: the
# equations and their solution, printed for comparison with PQCF (comparison/cumulants). Output: one JSON line.
using QuantumCumulants, ModelingToolkit, OrdinaryDiffEq
hc = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hc ⊗ ha
@qnumbers a::Destroy(h) σ::Transition(h)
@cnumbers Δ g κ γ ν
H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]
rates = [κ, γ, ν]
eqs = meanfield([a' * a, σ(:e, :e)], H, J; rates=rates, order=2)
eqs_c = complete(eqs)
println(stderr, eqs_c)
@named sys = ODESystem(eqs_c)
u0 = zeros(ComplexF64, length(eqs_c))
p0 = (Δ => 0.5, g => 1.5, κ => 1.0, γ => 0.25, ν => 4.0)
prob = ODEProblem(sys, u0, (0.0, 10.0), p0)
ts = [1.0, 2.0, 5.0, 10.0]
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, saveat=ts)
n = [real(sol[average(a' * a)][i]) for i in 1:length(ts)]
pe = [real(sol[average(σ(:e, :e))][i]) for i in 1:length(ts)]
println("{\"equations\": ", length(eqs_c), ", \"t\": ", ts, ", \"n\": ", n, ", \"pe\": ", pe, "}")
