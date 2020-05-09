"""
    setup_halofit(pk)

This function returns three parameterσ, [kσ, neff, C],  which are needed for computation of the halofit non-linear power spectrum.

*Reference*: Appendix C of Smith et al., MNRAS, 341, 1311 (2003)

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum **in units of (Mpc/h)^3** (NOTE h^-3!) with the argument k being the comoving wavenumber **in units of h/Mpc**.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function setup_halofit(pk)
    f(x) = √sigma2gaus(pk, x) - 1
    Rh = find_zero(f, 3)
    σ2 = sigma2gaus(pk, Rh)
    dlnσ2dlnRh = Rh * dsigma2gausdRh(pk, Rh) / σ2
    d2lnσ2dlnRh2 =
        dlnσ2dlnRh - dlnσ2dlnRh^2 + (Rh^2 / σ2) * d2sigma2gausdRh2(pk, Rh)
    kσ = 1/Rh
    neff = -3 - dlnσ2dlnRh
    C = -d2lnσ2dlnRh2
    return [kσ, neff, C]
end

"""
    halofit(pk, p, Ωmz, k)

This function returns a value of the non-linear power spectrum (**in units of Mpc^3/h^3**) at a specified value of the wavenumber **in units of h/Mpc**.

*Reference*: Appendix of Takahashi et al., ApJ, 761, 152 (2012)
- This reference updates parameters given in Appendix C of Smith et al., MNRAS, 341, 1311 (2003)

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum **in units of (Mpc/h)^3** (NOTE h^-3!) with the argument k being the comoving wavenumber **in units of h/Mpc**.
- `p::Array{T,1} where T <: Real`: an array containing [kσ, neff, C], which is obtained from ``p = setup_halofit(pk)``.
- `Ωmz::Real`: the matter density parameter at a given redshift z, i.e., Ωmz = Ωm(1+z)^3 / [Ωm(1+z)^3 + ΩΛ].
- `k::Real`: the comoving wavenumber **in units of h/Mpc**.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function halofit(pk, p::Array{T,1}, Ωmz::T, k::T) where T <: Real
    kσ, neff, C = p
    an =
        10^(
            1.5222 +
            2.8553neff +
            2.3706 * neff^2 +
            0.9903 * neff^3 +
            0.2250 * neff^4 - 0.6038C
        )
    bn = 10^(-0.5642 + 0.5864neff + 0.5716 * neff^2 - 1.5474C)
    cn = 10^(0.3698 + 2.0404neff + 0.8161 * neff^2 + 0.5869C)
    γn = 0.1971 - 0.0843neff + 0.8460C
    αn = abs(6.0835 + 1.3373neff - 0.1959 * neff^2 - 5.5274C)
    βn =
        2.0379 - 0.7354neff +
        0.3157 * neff^2 +
        1.2490 * neff^3 +
        0.3980 * neff^4 - 0.1682C
    μn = 0
    νn = 10^(5.2105 + 3.6902neff)
    f1 = Ωmz^-0.0307
    f2 = Ωmz^-0.0585
    f3 = Ωmz^0.0743
    y = k / kσ
    Δ2H = an * y^(3*f1) / (1 + bn * y^f2 + (cn * f3 * y)^(3 - γn))
    Δ2H /= 1 + μn / y + νn / y^2
    Δ2L = k^3 * pk.(k) / 2 / π^2
    Δ2Q = Δ2L * (1 + Δ2L)^βn / (1 + αn * Δ2L) * exp(-y / 4 - y^2 / 8)
    pk_nl = (Δ2Q + Δ2H) * 2 * π^2 / k^3
end
