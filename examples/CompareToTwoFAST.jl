using DoubleExponentialFormulas
using Plots, LaTeXStrings
using Dierckx
using DelimitedFiles, Printf
using TwoFAST

# %% Specify a redshift
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
ired = 2
z = redshift[ired]

# %% Read in linear power spectra computed by CLASS and spline interpolate in k
filename = @sprintf("data/matterpower05_z%1d_pk.dat", ired)
d = readdlm(filename, comments = true)
pk = Spline1D(d[:, 1], d[:, 2])

# %% Compute correlation function with `quaddeo`
kmin, kmax = 1e-4, Inf
R1 = Float64.(10:2:200)
ξ1 = zeros(length(R1))
for i = 1:length(R1)
   func1(k) = k^2 * pk(k) * sin(k * R1[i]) / (k * R1[i]) / 2 / π^2
   ξ1[i], err = quaddeo(func1, R1[i], 0, kmin, kmax)
end

# %% Compute correlation function with `xicalc` in TwoFAST
R2, ξ2 = xicalc(pk, 0, 0)

# %% Plot results and save to "compare_xi_twofast.pdf"
p = plot(
   R1,
   R1 .^ 2 .* ξ1,
   m = 2,
   c = :red,
   lab = "Computed with quaddeo",
   ylab = L"R^2\xi(R)",
   xlab = L"R~[Mpc/h]",
   title = "Redshift = $z",
   xlims = [0,200],
)
p = plot!(R2, R2 .^ 2 .* ξ2, ls = :dash, c = :red, lab = "with TwoFAST")
savefig("compare_xi_twofast.pdf")
display(p)
