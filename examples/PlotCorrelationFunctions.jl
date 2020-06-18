using MatterPower
using DoubleExponentialFormulas
using Plots, LaTeXStrings
using Dierckx
using DelimitedFiles, Printf

# %% Specify a redshift
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
ired = 1

# %% Read in linear and non-linear power spectra computed by CLASS and spline interpolate in k
filename = @sprintf("data/matterpower05_z%1d_pk.dat", ired)
d = readdlm(filename, comments = true)
pklin_class = Spline1D(d[:, 1], d[:, 2])
filename = @sprintf("data/matterpower05_z%1d_pk_nl.dat", ired)
d = readdlm(filename, comments = true)
pknl_class = Spline1D(d[:, 1], d[:, 2])

# %% Define a function to return a linear matter power spectrum (in units of Mpc^3/h^3)
# as a function of the comoving wavenumber, kovh, in units of h/Mpc.
# Here is an example using Einstein & Hu's analytical transfer function

# Cosmological parameters
As, ns, kpivot = 2.097e-9, 0.9652, 0.05
Ωm, ΩΛ, Ωb = 0.315, 0.685, 0.049
Ωk = 1 - Ωm - ΩΛ
h0 = 0.674
ωm, fb = Ωm * h0^2, Ωb / Ωm

# Tabulate linear growth factor as a function of scale factor, a
sol = setup_growth(Ωm, ΩΛ)
a = 1 / (1 + redshift[ired])
D1 = sol(a)[1]

pk(kovh) =
   D1^2 *
   As *
   (kovh * h0 / kpivot)^(ns - 1) *
   (2 * kovh^2 * 2998^2 / 5 / Ωm)^2 *
   t_nowiggle(kovh * h0, ωm, fb)^2 *
   2 *
   π^2 / kovh^3

# %% Alternatively you may read in pre-computed data and define a spline function
# using Dierckx
# pk = Spline1D(tabulated_kovh, tabulated_power_spectrum)

# %% Compute the non-linear power spectrum using halofit
# Get three parameters [kσ, neff, C] needed for the halofit transform
p = kσ, neff, C = setup_halofit(pk)

# Calculate Ωm(z) = Ωm(z=0)(1+z)^3 / E^2(z), where E^2(z) = H^2(z) / H0^2
z = redshift[ired]
E2 = Ωm * (1 + z)^3 + Ωk * (1 + z)^2 + ΩΛ
Ωmz = Ωm * (1 + z)^3 / E2

# Define a function to return a halofit non-linear power spectrum
pknl(kovh) = halofit(pk, p, Ωmz, kovh) # Mpc^3/h^3

# %% Compute correlation functions
kmin, kmax = 1e-4, Inf
R = Float64.(10:2:200)
ξ1 = zeros(length(R))
ξ2 = zeros(length(R))
ξ3 = zeros(length(R))
ξ4 = zeros(length(R))
for i = 1:length(R)
   func1(k) = k^2 * pk(k) * sin(k * R[i]) / (k * R[i]) / 2 / π^2
   ξ1[i], err = quaddeo(func1, R[i], 0, kmin, kmax)
   func2(k) = k^2 * pknl(k) * sin(k * R[i]) / (k * R[i]) / 2 / π^2
   ξ2[i], err = quaddeo(func2, R[i], 0, kmin, kmax)
   func3(k) = k^2 * pklin_class(k) * sin(k * R[i]) / (k * R[i]) / 2 / π^2
   ξ3[i], err = quaddeo(func3, R[i], 0, kmin, kmax)
   func4(k) = k^2 * pknl_class(k) * sin(k * R[i]) / (k * R[i]) / 2 / π^2
   ξ4[i], err = quaddeo(func4, R[i], 0, kmin, kmax)
end

# %% Plot results and save to "xi.pdf"
p = plot(
   R,
   R .^ 2 .* ξ4,
   m = 2,
   c = :red,
   lab = "Halofit ξ(R) with BAO",
   ylab = L"R^2\xi(R)",
   xlab = L"R~[Mpc/h]",
   title = "Redshift = $z",
)
p =
   plot!(R, R .^ 2 .* ξ2, ls = :dash, c = :red, lab = "Halofit ξ(R) w/o BAO")
p = plot!(R, R .^ 2 .* ξ3, m = 2, c = 1, lab = "Linear ξ(R) with BAO")
p = plot!(R, R .^ 2 .* ξ1, ls = :dash, c = 1, lab = "Linear ξ(R) w/o BAO")
savefig("xi.pdf")
display(p)
