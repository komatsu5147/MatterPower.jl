using MatterPower
using PyCall
using Roots
using Tables, CSV
# %% Call the python wrapper for CLASS, `classy`, via PyCall
classy = pyimport("classy")
# Create an instance of the CLASS wrapper
cosmo = classy.Class()
# Create a dictionary of the cosmological parameters
deg_ncdm = true # 3 massive neutrinos with degenerate mass?
params = Dict(
   "output" => "mPk",
   "P_k_max_h/Mpc" => 100,
   "z_pk" => 5,
   "non linear" => "halofit",
   "A_s" => 2.097e-9,
   "n_s" => 0.9652,
   "k_pivot" => 0.05,
   "h" => 0.6737,
   "omega_b" => 0.02233,
   "omega_cdm" => 0.1198,
   "N_ncdm" => 1,
)
if deg_ncdm
   push!(params, "m_ncdm" => 0.02, "deg_ncdm" => 3, "N_ur" => 0.00641)
   mν = params["m_ncdm"] * params["deg_ncdm"]
else
   push!(params, "m_ncdm" => 0.06, "N_ur" => 2.0328)
   mν = params["m_ncdm"]
end
h0 = params["h"]
Ωcb = (params["omega_b"] + params["omega_cdm"]) / h0^2
# Set the parameters and run the CLASS code
cosmo.set(params)
cosmo.compute()

# %% Compute Ωgrav and Ωtherm at seven redshifts
redshift = 0:0.1:2
nred = length(redshift)
Rstar = zeros(nred)
Mstar = zeros(nred)
for ired = 1:nred
   z = redshift[ired]
   # Define functions to return a linear baryon+CDM power spectrum
   # Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
   # return power spectra in units of Mpc^3 (no 1/h^3).
   pkcb(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
   # %% Compute the non-linear mass, M*, defined by σ(R*,z) = 1.6865
   # and M* = (4π/3)ρm R*^3 = (4π/3)ρc Ωm R*^3
   # ρc is the present-day critical density of the Universe in units of h^2 M⊙/Mpc^3
   f(R) = √sigma2(pkcb, R) - 1.6865
   Rstar[ired] = find_zero(f, (0.2, 3), Bisection())
end
ρc = 2.775e11
Mstar = 4π / 3 * ρc * Ωcb * Rstar .^ 3

#%% Save table data to table1.csv
t = Tables.table([redshift Rstar Mstar])
header = ["z", "Rstar", "Mstar"]
CSV.write("redshift_Rstar_Mstar.csv", t, header = header)

# %% Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
cosmo.struct_cleanup()
# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not neededc
cosmo.empty()
