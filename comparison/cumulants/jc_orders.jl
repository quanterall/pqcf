# The driven, damped Jaynes-Cummings model at cumulant orders 2-4 in QuantumCumulants.jl: <ad a> at t = 2 from vacuum,
# for comparison with PQCF (comparison/cumulants). Output: one JSON line.
using QuantumCumulants, ModelingToolkit, OrdinaryDiffEq
hc = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hc ⊗ ha
@qnumbers a::Destroy(h) σ::Transition(h)
@cnumbers g E κ γ
H = g * (a' * σ(:g, :e) + a * σ(:e, :g)) + E * (a + a')
J = [a, σ(:g, :e)]
rates = [κ, γ]
out = Dict{Int,Any}()
for order in 2:4
    eqs = complete(meanfield([a' * a], H, J; rates=rates, order=order))
    @named sys = ODESystem(eqs)
    u0 = zeros(ComplexF64, length(eqs))
    prob = ODEProblem(sys, u0, (0.0, 2.0), (g => 1.0, E => 0.5, κ => 0.5, γ => 0.25))
    sol = solve(prob, Vern9(); abstol=1e-14, reltol=1e-14, saveat=[2.0])
    out[order] = (length(eqs), real(sol[average(a' * a)][end]))
end
println("{", join(["\"$k\": {\"equations\": $(out[k][1]), \"n\": $(out[k][2])}" for k in 2:4], ", "), "}")
