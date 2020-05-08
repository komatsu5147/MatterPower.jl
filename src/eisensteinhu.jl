"""
    t_nowiggle(k, ωm, fbaryon)

Transfer function for the linear matter power spectrum without Baryon Acoustic Oscillation (the so-called *no-wiggle* transfer function).

*Reference*: Equation (29-31) of Eisenstein and Hu, ApJ, 496, 605 (1998)

# Arguments
- `k`: the comoving wavenumber **in units of 1/Mpc**. I.e., there is no *h* in the wavenumber.
- `ωm`: the physical baryon density parameter, `ωm` = Ωm h^2
- `fbaryon`: the baryon fractionr,  `fbaryon` = Ωb/Ωb

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function t_nowiggle(k, ωm, fb) # k is in units of Mpc^{-1} [without h]
    α = 1 - 0.328fb * log(431ωm) + 0.38 * fb^2 * log(22.3ωm)
    s = 44.5 * log(9.83 / ωm) / √(1 + 10 * (fb * ωm)^(3 / 4))
    Γ = ωm * (α + (1 - α) / (1 + (0.43 * k * s)^4))
    q = k * (2.725 / 2.7)^2 / Γ
    log(2 * exp(1) + 1.8q) /
    (log(2 * exp(1) + 1.8q) + (14.2 + 731 / (1 + 62.5q)) * q^2)
end
