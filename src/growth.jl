"""
    setup_growth(Ωm, ΩΛ[, a1=1e-2])

Return the solutions of coupled differential equations given in `growth!`,
with the initial condition such that ``δ(a) = a`` at `a = a1`. The default value is `a1 = 1e-2`.

# Arguments
- `Ωm::Real`: present-day total matter density parameter.
- `ΩΛ::Real`: present-day dark energy density parameter (for cosmological constant).

# Optional arguments
- `a1::Real=1e-2`: initial scale factor, at which the initial condition is set as δ(a1) = a1.

# Examples
```julia-repl
julia> Ωm, ΩΛ = 0.3, 0.7
(0.3, 0.7)

julia> sol = MatterPower.setup_growth(Ωm, ΩΛ)
retcode: Success
Interpolation: specialized 4th order "free" interpolation
t: 14-element Array{Float64,1}:
 0.01
 0.022952118470240847
 0.03152976635130906
 0.04919933323749941
 ⋮
 0.7292247303730331
 0.994509796947469
 1.0
u: 14-element Array{Array{Float64,1},1}:
 [0.01, -0.05477231965144438]
 [0.022953661875878814, -0.08298398087457971]
 [0.0315315317884963, -0.09726149943169236]
 [0.049200231832575124, -0.12149028570497]
 ⋮
 [0.6428440417117532, -0.39916939875846624]
 [0.7768316486186253, -0.3997506971222691]
 [0.7790350080966789, -0.3994861471523371]

julia> redshift = 1.2
1.2

julia> a = 1/(1+redshift)
0.45454545454545453

julia> δ = sol(a)[1]
0.43811688788577136

julia> θ = sol(a)[2]
-0.3526347043144974
```
"""
function setup_growth(Ωm::Real, ΩΛ::Real, a1::Real=1e-2)
    Ωk = 1 - Ωm - ΩΛ
    E1 = √(Ωm / a1^3 + Ωk / a1^2 + ΩΛ)
    u0 = [a1; -a1^2 * E1]
    tspan = (a1, 1.0)
    prob = ODEProblem(growth!, u0, tspan, [Ωm, ΩΛ])
    sol = solve(prob, Tsit5())
end

"""
    growth!(du, u, p, a)

Coupled differential equations to obtain the growth factor of linear density fluctuations.

1. Continuity Equation, ``δ' + θ = 0``

2. Euler equation, ``θ' + (a'/a)θ - k^2 ϕ = 0``
3. Poisson equation, ``k^2ϕ = -4πG a^2 ρm δ``

where θ = div(velocity), ϕ is the Newtonian gravitational potential, and the primes denote conformal time derivatives.

To solve these equations, we rewrite them using derivatives with respect to the scale factor, `a`.
1. `du[1]` = ``dδ/da = -θ/a^2/E(a)``
2. `du[2]` = ``dθ/da = -θ/a + q^2ϕ/a^2/E(a)``
3. ``q^2ϕ = -(3/2)Ωm δ/a``

where

- `u[1]` = δ, `u[2]` = θ = div(velocity)/H0
- `p[1]` = Ωm, `p[2]` = ΩΛ
- ``E(a) = H(a)/H0 = √(Ωm/a^3 + Ωk/a^2 + ΩΛ)``
- ``q = k/H0``

These equations will be solved by Julia's ODE solver as

`ODEProblem(growth!, u0, tspan, p)`

See this [documentation](https://github.com/SciML/OrdinaryDiffEq.jl) for how to use this ODE solver.
"""
function growth!(du, u, p, a)
    δ, θ = u
    Ωm, ΩΛ = p
    Ωk = 1 - Ωm - ΩΛ
    E = √(Ωm / a^3 + Ωk / a^2 + ΩΛ) # E(a) = H(a)/H0
    q2ϕ = -(3 / 2) * Ωm * δ / a # Poisson equation
    du[1] = dδda = -θ / a^2 / E # Continuity equation
    du[2] = dθda = -θ / a + q2ϕ / a^2 / E # Euler equation
end
