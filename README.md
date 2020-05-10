# MatterPower

This package contains functions to compute
1. Growth factor of the linear matter density fluctuation for ΛCDM cosmology
2. "No-wiggle" transfer function for the linear matter power spectrum
3. Non-linear matter power spectrum
4. Variance of the matter density fluctuation and its derivatives

The package contains
1. Linear growth factor
  - `setup_growth(Ωm, ΩΛ)`: returns an interpolation function containing solutions of the linear matter density (δ) and velocity divergence (θ/H0) fields as a function of the scale factor, a  
2. Linear matter transfer function
  - `t_nowiggle(k, ωm, fbaryon)`: returns a "no-wiggle transfer function", T0(k), from Equation (29-31) of [Eisenstein and Hu, ApJ, 496, 605 (1998)](https://iopscience.iop.org/article/10.1086/305424)
3. Non-linear matter power spectrum
  - `setup_halofit(pk)`: returns an array containing three variables `p=[kσ, neff, C]` defined in Equation (C6-C8) of [Smith et al., MNRAS, 341, 1311 (2003)](https://academic.oup.com/mnras/article/341/4/1311/1038529)
  - `halofit(pk, p, Ωmz, k)`: returns a "halofit" non-linear matter power spectrum at a specified value of the comoving wavenuember `k` using [Smith et al.](https://academic.oup.com/mnras/article/341/4/1311/1038529)'s method with the parameters given in Appendix of [Takahashi et al., ApJ, 761, 152 (2012)](https://iopscience.iop.org/article/10.1088/0004-637X/761/2/152)
    - The terms proportional to (1+w) in Equation (A6, A7) of Takahashi et al. are not included in this function
4. Variance of the matter density fluctuation
  - `sigma2(pk, R)`: returns variance of the matter density fluctuation, σ^2, smoothed by a top-hat filter with radius `R`
  - `dsigma2dR(pk, R)`: returns a derivative of variance, dσ^2/dR, with respect to a top-hat filter with radius `R`
  - `sigma2gaus(pk, R)`: returns variance of the matter density fluctuation, σ^2, smoothed by a Gaussian filter with width `R`
  - `dsigma2gausdR(pk, R)`: returns a derivative of variance, dσ^2/dR, with respect to a Gaussian filter with width `R`
  - `d2sigma2gausdR2(pk, R)`: returns a secon derivative of variance, d^2σ^2/dR^2, with respect to a Gaussian filter with width `R`

## Arguments
- `Ωm::Real`: present-day total matter density parameter.
- `ΩΛ::Real`: present-day dark energy density parameter (for cosmological constant).
- `k::Real`: comoving wavenumber.
- `ωm::Real`: physical baryon density parameter, ``ωm = Ωm h^2``
- `fbaryon::Real`: baryon fraction,  ``fbaryon = Ωb/Ωb``.
- `pk::Any`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber. This can be an interpolation function constructed from tabulated data.
- `p::Array{T,1} where T <: Real`: an array containing three variables [kσ, neff, C], which are obtained from `p = setup_halofit(pk)`.
- `Ωmz::Real`: matter density parameter at a given redshift z, e.g., ``Ωmz = Ωm(1+z)^3 / [Ωm(1+z)^3 + Ωk(1+z)^2 + ΩΛ]``.
- `R::Real`: top-hat or Gaussian smoothing scale.

## Motivation

The linear matter power spectrum computed from the "no-wiggle" transfer function `t_nowiggle(k, ωm, fbaryon)` is not very accurate because it does not contain the Baryon Acoustic Oscillation (BAO) and it is a fitting function. If you need a very accurate linear matter power spectrum, you should use linear Boltzmann solvers such as [CAMB](https://github.com/cmbant/CAMB) and [CLASS](https://github.com/lesgourg/class_public).

Nevertheless, if you do not need BAO, this linear matter power spectrum is reasonably accurate: it achieves precision better than 1.5% at all wavenumbers below k < 10 h/Mpc for ΛCDM cosmology with the standard thermal history of the Universe. This is quite sufficient for many purposes in the cosmological research (unless you need BAO, of course). Therefore, you can calculate many quantities in your research with decent accuracy very quickly using the functions given in this package.

Finding new, (hopefully) interesting problems is the hardest part of research. To explore new ideas, precision is often not needed; thus, you can use the functions given in this package to explore new ideas. Once you find interesting results that are worth exploring further with better precision, you can update your calculations with CAMB or CLASS.

Also, having a single, self-contained, simple and well-documented package should be useful for educational purposes. My hope is that beginning students would take a look at this package and gain insights into how the growth of linear density fluctuation and linear and non-linear matter power spectra are computed, before using public linear Boltzmann solvers as a black box.

Having said this, the non-linear matter power spectrum functions `setup_halofit(pk)` and `halofit(pk, p, Ωmz, k)` in this package can take any linear power spectrum inputs `pk` you wish; i.e, you do not have to use `t_nowiggle(k, ωm, fbaryon)`, but can supply linear power spectra computed from CAMB and CLASS.

## How the linear growth factor is calculated

Growth of the linear matter density fluctuation, δ, is determined by three equations: (1) continuity, (2) Euler, and (3) Poisson equations.

1. Continuity Equation, ``δ' + θ = 0``
2. Euler equation, ``θ' + (a'/a)θ - k^2 ϕ = 0``
3. Poisson equation, ``k^2ϕ = -4πG a^2 ρm δ``

where θ = div(velocity), ϕ is the Newtonian gravitational potential, and the primes denote conformal time derivatives.

To solve these equations, we rewrite them using derivatives with respect to the scale factor, `a`.

1. ``dδ/da = -θ/a^2/E(a)``
2. ``dθ/da = -θ/a + q^2ϕ/a^2/E(a)``
3. ``q^2ϕ = -(3/2)Ωm δ/a``

where
- ``θ = div(velocity)/H0``
- ``E(a) = H(a)/H0 = √(Ωm/a^3 + Ωk/a^2 + ΩΛ)``
- ``q = k/H0``

These equations are put in the form of a function `growth!(du, u, p, a)`
```
function growth!(du, u, p, a)
    δ, θ = u
    Ωm, ΩΛ = p
    Ωk = 1 - Ωm - ΩΛ
    E = √(Ωm / a^3 + Ωk / a^2 + ΩΛ) # E(a) = H(a)/H0
    k2ϕ = -(3 / 2) * Ωm * δ / a # Poisson equation
    du[1] = dδda = -θ / a^2 / E # Continuity equation
    du[2] = dθda = -θ / a + k2ϕ / a^2 / E # Euler equation
end
```
and solved by Julia's ODE solver as
```
prob = ODEProblem(growth!, u0, tspan, [Ωm, ΩΛ])
sol = solve(prob, Tsit5())
```
with
- `du[1]` = dδ/da, `du[2]` = dθ/da
- `u[1]` = δ, `u[2]` = θ
- `p[1]` = Ωm, `p[2]` = ΩΛ
- `tspan = (a1, 1.0)` with `a1 = 0.01`
- `u0 = [a1; -a1^2 * E(a1)]`

See this [documentation](https://github.com/SciML/OrdinaryDiffEq.jl) for how to use this ODE solver.

## Example Juia code

The following example code (avaiable in [examples/PowerSpectrum.jl](https://github.com/komatsu5147/MatterPower.jl/blob/master/examples/PowerSpectrum.jl)) computes the linear growth factor, the linear and non-linear power spectra, the r.m.s. mass density fluctuation σ(R), and the non-linear mass M* defined by σ(M*) = 1.6865.

If you would like to generate a nice figure showing linear and non-linear P(k) as a function of k, take a look at [examples/PlotPowerSpectrum.jl](https://github.com/komatsu5147/MatterPower.jl/blob/master/examples/PowerSpectrum.jl).
```
using MatterPower
using Roots

# %% Specify a redshift
redshift = 0

# %% Define a function to return a linear matter power spectrum (in units of Mpc^3/h^3)
# as a function of the comoving wavenumber, k_ov_h, in units of h/Mpc.
# Here is an example using Einstein & Hu's analytical transfer function

# Cosmological parameters
As, ns, kpivot = 2.097e-9, 0.9652, 0.05
Ωm, ΩΛ, Ωb = 0.315, 0.685, 0.049
Ωk = 1 - Ωm - ΩΛ
h0 = 0.674
ωm, fb = Ωm * h0^2, Ωb / Ωm

# Tabulate linear growth factor as a function of scale factor, a
sol = setup_growth(Ωm, ΩΛ)
a = 1 / (1 + redshift)
D1 = sol(a)[1]

pk(k_ov_h) =
   D1^2 *
   As *
   (k_ov_h * h0 / kpivot)^(ns - 1) *
   (2 * k_ov_h^2 * 2998^2 / 5 / Ωm)^2 *
   t_nowiggle(k_ov_h * h0, ωm, fb)^2 *
   2 *
   π^2 / k_ov_h^3

# %% Alternatively you may read in pre-computed data and define a spline function
# using Dierckx
# pk = Spline1D(tabulated_k_ov_h, tabulated_power_spectrum)

# %% Compute the r.m.s. mass fluctuation with a top-hat radius Rh
Rh = 8 # Mpc/h
σ2 = sigma2(pk, Rh)
println("σ(R) = ", √σ2, " for R = ", Rh, " Mpc/h")

# %% Compute the non-linear mass, M*, defined by σ(R*,z) = 1.6865
# and M* = (4π/3)ρm R*^3 = (4π/3)ρc Ωm R*^3
# ρc is the present-day critical density of the Universe in units of h^2 M⊙/Mpc^3
ρc = 2.775e11
f(x) = √sigma2(pk, x) - 1.6865
Rstar = find_zero(f, (1e-2, 1e2), Bisection())
println("R* = ", Rstar, " h/Mpc")
Mstar = 4π / 3 * ρc * Ωm * Rstar^3
println("Non-linear mass M* = ", Mstar, " M⊙/h")

# %% Compute the non-linear power spectrum using halofit
# Get three parameters [kσ, neff, C] needed for the halofit transform
p = kσ, neff, C = setup_halofit(pk)
# Calculate Ωm(z) = Ωm(z=0)(1+z)^3 / E^2(z), where E^2(z) = H^2(z) / H0^2
z = redshift
E2 = Ωm * (1 + z)^3 + Ωk * (1 + z)^2 + ΩΛ
Ωmz = Ωm * (1 + z)^3 / E2
# Get a halofit non-linear power spectrum at a specified value of the comoving wavenumber
k_ov_h = 1 # h/Mpc
pknl = halofit(pk, p, Ωmz, k_ov_h) # Mpc^3/h^3
println("k = ", k_ov_h, " h/Mpc:")
println("Non-linear P(k) = ", pknl, " Mpc^3/h^3")
println("Linear P(k) = ", pk(k_ov_h), " Mpc^3/h^3")
```
## Acknowledgment

The functions provided in this package are adapted from [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/), which is based upon work supported in part by NSF under Grant AST-0807649 and PHY-0758153, NASA under Grant NNX08AL43G, and Alfred P. Sloan Research Foundation via a Sloan Fellowship. This work is also supported in part by JSPS KAKENHI Grant Number JP15H05896.
