using MatterPower
using PyCall
using Plots, LaTeXStrings

# %% Specify a redshift
redshift = 0

# %% Define functions to return linear and non-linear matter power spectra (in units of Mpc^3/h^3)
classy = pyimport("classy")
cosmo = classy.Class()

# Set the parameters to the CLASS code
params = Dict(
    "output" => "mPk",
    "P_k_max_h/Mpc" => 30,
    "non linear" => "halofit",
    "A_s" => 2.097e-9,
    "n_s" => 0.9652,
    "k_pivot" => 0.05,
    "h" => 0.674,
    "Omega_b" => 0.049,
    "Omega_cdm" => 0.266,
)

# Set the parameters and run the CLASS code
cosmo.set(params)
cosmo.compute()

# Define functions to return linear and non-linear power spectra
# Note: CLASS codes take wavenumbers in units of 1/Mpc (no h) and
# return power spectra in units of Mpc^3 (no 1/h^3).
pklin_class(kovh) = cosmo.pk_lin(kovh * h0, redshift) * h0^3
pknl_class(kovh) = cosmo.pk(kovh * h0, redshift) * h0^3

# %% Plot results and save to "classpk.pdf"
lnk = log(3e-4):0.3:log(30)
kovh = exp.(lnk)
p = plot(
    kovh,
    pknl_class.(kovh),
    xaxis = :log,
    yaxis = :log,
    m = 2,
    lab = "Non-linear P(k)",
    ylab = L"P(k)~[h^{-3}~Mpc^3]",
    xlab = L"k~[h~Mpc^{-1}]",
)
p = plot!(kovh, pklin_class.(kovh), ls = :dash, lab = "Linear P(k)")
savefig("classpk.pdf")
display(p)
