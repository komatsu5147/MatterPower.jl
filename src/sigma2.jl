"""
    sigma2(pk, R)

Variance of mass fluctuation, σ^2, with a top-hat filter with radius `R`

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `R::Real`: a top-hat smoothing scale.
   - **R times k must be dimensionless**. For example, if k is in units of h/Mpc, `R` must be in units of Mpc/h.

`sigma2` is computed as

``σ^2(R) = ∫_0^∞ k^2dk W^2(kR) pk(k) / (2π^2)``

where and W(kR) is a window function given by

``W(x) = (3/x)j_1(x) = (3/x)(sin(x)/x^2 - cos(x)/x)``

j_1(x) is the spherical Bessel function of order 1.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function sigma2(pk, R::Real)
   function dσ2dk(k)
      x = k * R
      W = (3 / x) * (sin(x) / x^2 - cos(x) / x)
      dσ2dk = W^2 * pk(k) * k^2 / 2 / π^2
   end
   res, err = quadde(dσ2dk, 0, Inf)
   σ2 = res
end

"""
    dsigma2dR(pk, R)

Derivative of variance of mass fluctuation, dσ^2/dR,  with respect to a top-hat filter with radius `R`.

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `R::Real`: a top-hat smoothing scale.
   - **R times k must be dimensionless**. For example, if k is in units of h/Mpc, `R` must be in units of Mpc/h.

`dsigma2/dR` is computed as

``dσ^2/dR = 2 ∫_0^∞ k^2dk W(kR) dW/dR pk(k) / (2π^2)``

where and W(kR) is a window function given by

``W(x) = (3/x)j_1(x) = (3/x)(sin(x)/x^2 - cos(x)/x)``

and dW/dR is its derivative:

``dW/dR = k dW/dx = (-3k/x)j_2(x) = (-3k/x)[(3/x^2-1)sin(x)/x - 3cos(x)/x^2]``

j_n(x) is the spherical Bessel function of order n.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function dsigma2dR(pk, R::Real)
   function dσ2dRdk(k)
      x = k * R
      W = (3 / x) * (sin(x) / x^2 - cos(x) / x)
      dWdx = (-3 / x) * ((3 / x^2 - 1) * sin(x) / x - 3 * cos(x) / x^2)
      dσ2dRdk = 2 * W * dWdx * pk(k) * k^3 / 2 / π^2
   end
   res, err = quadde(dσ2dRdk, 0, Inf)
   dσ2dR = res
end

"""
    sigma2gaus(pk, R)

Variance of mass fluctuation, σ^2(`R`), with a Gaussian filter with width `R`.

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `R::Real`: a top-hat smoothing scale.
   - **R times k must be dimensionless**. For example, if k is in units of h/Mpc, `R` must be in units of Mpc/h.

`sigma2gaus(R)` is computed as

``σ^2(R) = ∫_0^∞ k^2dk exp(-k^2 R^2) pk(k) / (2π^2)``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function sigma2gaus(pk, R::Real)
   function dσ2dk(k)
      x = k * R
      dσ2dk = exp(-x^2) * pk(k) * k^2 / 2 / π^2
   end
   res, err = quadde(dσ2dk, 0, Inf)
   σ2 = res
end

"""
    dsigma2gausdR(pk, R)

Derivative of variance of mass fluctuation, dσ^2/dR, with respect to a Gaussian filter with width `R`.

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `R::Real`: a top-hat smoothing scale.
   - **R times k must be dimensionless**. For example, if k is in units of h/Mpc, `R` must be in units of Mpc/h.

`dsigma2gaus/dR` is computed as

``dσ^2/dR = -2R ∫_0^∞ k^4dk exp(-k^2 R^2) pk(k) / (2π^2)``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function dsigma2gausdR(pk, R::Real)
   function dσ2dRdk(k)
      x = k * R
      dσ2dRdk = -2R * exp(-x^2) * pk(k) * k^4 / 2 / π^2
   end
   res, err = quadde(dσ2dRdk, 0, Inf)
   dσ2dR = res
end

"""
    d2sigma2gausdR2(pk, R)

Second derivative of variance of mass fluctuation, d^2σ^2/dR^2, with respect to a Gaussian filter with width `R`.

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `R::Real`: a top-hat smoothing scale.
   - **R times k must be dimensionless**. For example, if k is in units of h/Mpc, `R` must be in units of Mpc/h.

`d2sigma2gaus/dR2` is computed as

``d^2σ^2/dR^2 = -2 ∫_0^∞ (1 - 2x^2) k^4dk exp(-k^2 R^2) pk(k) / (2π^2)``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function d2sigma2gausdR2(pk, R::Real)
   function d2σ2dR2dk(k)
      x = k * R
      d2σ2dR2dk = (-2 + 4 * x^2) * exp(-x^2) * pk(k) * k^4 / 2 / π^2
   end
   res, err = quadde(d2σ2dR2dk, 0, Inf)
   d2σ2dR2 = res
end
