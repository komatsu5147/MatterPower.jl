using MatterPower
using Roots
using Plots, LaTeXStrings
using Dierckx
using DelimitedFiles, Printf

# %% Loop over redshifts
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
for ired = 1:7
   # %% Read in linear and non-linear power spectra computed by CLASS
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

   # %% Compute the r.m.s. mass fluctuation with a top-hat radius Rh
   Rh = 8 # Mpc/h
   σ2 = sigma2(pk, Rh)
   println("σ(R) = ", √σ2, " for R = ", Rh, " Mpc/h")

   # %% Compute the non-linear mass, M*, defined by σ(R*,z) = 1.6865
   # and M* = (4π/3)ρm R*^3 = (4π/3)ρc Ωm R*^3
   # ρc is the present-day critical density of the Universe in units of h^2 M⊙/Mpc^3
   ρc = 2.775e11
   f(x) = √sigma2(pk, x) - 1.6865
   Rstar = find_zero(f, (1e-2, 3), Bisection())
   println("R* = ", Rstar, " h/Mpc")
   Mstar = 4π / 3 * ρc * Ωm * Rstar^3
   println("Non-linear mass M* = ", Mstar, " M⊙/h")

   # %% Compute the non-linear power spectrum using halofit
   # Get three parameters [kσ, neff, C] needed for the halofit transform
   p = kσ, neff, C = setup_halofit(pk)

   # Calculate Ωm(z) = Ωm(z=0)(1+z)^3 / E^2(z), where E^2(z) = H^2(z) / H0^2
   z = redshift[ired]
   E2 = Ωm * (1 + z)^3 + Ωk * (1 + z)^2 + ΩΛ
   Ωmz = Ωm * (1 + z)^3 / E2

   # Define a function to return a halofit non-linear power spectrum
   pknl(kovh) = halofit(pk, p, Ωmz, kovh) # Mpc^3/h^3

   # %% Plot results and save to "comparison_z[1-7].pdf"
   lnk = log(3e-4):1e-2:log(30)
   kovh = exp.(lnk)
   p = plot(
      kovh,
      pknl.(kovh) ./ pknl_class.(kovh) .- 1,
      xaxis = :log,
      lab = "Non-linear P(k)",
      ylab = L"P_{MatterPower}(k)/P_{CLASS}(k) - 1",
      xlab = L"k~[h~Mpc^{-1}]",
      title = @sprintf("Redshift: z = %1.1f", redshift[ired]),
   )
   p = plot!(
      kovh,
      pk.(kovh) ./ pklin_class.(kovh) .- 1,
      ls = :dash,
      lab = "Linear P(k)",
   )
   filename = @sprintf("comparison_z%1d.pdf", ired)
   savefig(filename)
   display(p)

   # %% Compute the halofit power spectrum from CLASS's linear power spectrum
   p = kσ, neff, C = setup_halofit(pklin_class)
   pknl2(kovh) = halofit(pklin_class, p, Ωmz, kovh) # Mpc^3/h^3

   # %% Plot results and save to "halofit_comparison_z[1-7].pdf"
   p = plot(
      kovh,
      pknl2.(kovh) ./ pknl_class.(kovh) .- 1,
      xaxis = :log,
      lab = "Non-linear P(k) from CLASS's linear P(k) and halofit",
      ylab = L"P_{MatterPower}(k)/P_{CLASS}(k) - 1",
      xlab = L"k~[h~Mpc^{-1}]",
      title = @sprintf("Redshift: z = %1.1f", redshift[ired]),
      legend = :bottomleft,
   )
   filename = @sprintf("halofit_comparison_z%1d.pdf", ired)
   savefig(filename)
   display(p)
end
