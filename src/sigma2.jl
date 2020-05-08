"""
    sigma2(pk, Rh)

Variance of mass fluctuation, sigma^2(`Rh`).

# Arguments
- `pk`(k): a function which returns the linear matter power spectrum **in units of (Mpc/h)^3** (NOTE h^-3!) with the argument k being the comoving wavenumber **in units of h/Mpc**.
- `Rh`: a top-hat smoothing scale **in units of Mpc/h** (NOTE h^-1!).

`sigma2(Rh)` is computed as

``sigma^2(Rh) = ∫_0^∞ k^2dk W^2(kRh) pk(k) / (2π^2)``

where and W(kRh) is a window function given by

``W(x) = (3/x)j_1(x) = (3/x)(sin(x)/x^2 - cos(x)/x)``

j_1(x) is the spherical Bessel function of order 1.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function sigma2(pk, Rh)
   function dσ2dk(k)
      x = k * Rh  # Rh is in units of Mpc/h
      W = (3 / x) * (sin(x) / x^2 - cos(x) / x)
      dσ2dk = W^2 * pk(k) * k^2 / 2 / π^2 # pk is in units of Mpc^3/h^3
   end
   res, err = hquadrature(dσ2dk, 0, 20 / Rh)
   return res
end

"""
    dsigma2dRh(pk, Rh)

Derivative of variance of mass fluctuation with respect to `Rh`, dsigma^2/dRh.

# Arguments
- `pk`(k): a function which returns the linear matter power spectrum **in units of (Mpc/h)^3** (NOTE h^-3!) with the argument k being the comoving wavenumber **in units of h/Mpc**.
- `Rh`: a top-hat smoothing scale **in units of Mpc/h** (NOTE h^-1!).

`dsigma^2/dRh` is computed as

``dsigma^2/dRh = 2 ∫_0^∞ k^2dk W(kRh) dW/dRh pk(k) / (2π^2)``

where and W(kRh) is a window function given by

``W(x) = (3/x)j_1(x) = (3/x)(sin(x)/x^2 - cos(x)/x)``

and dW/dRh is its derivative:

``dW/dRh = k dW/dx = (-3k/x)j_2(x) = (-3k/x)[(3/x^2-1)sin(x)/x - 3cos(x)/x^2]``

j_n(x) is the spherical Bessel function of order n.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function dsigma2dRh(pk, Rh)
   function dσ2dRdk(k) # k is in units of h/Mpc
      x = k * Rh  # Rh is in units of Mpc/h
      W = (3 / x) * (sin(x) / x^2 - cos(x) / x)
      dWdx = (-3 / x) * ((3 / x^2 - 1) * sin(x) / x - 3 * cos(x) / x^2)
      dσ2dRdk = 2 * W * dWdx * pk(k) * k^3 / 2 / π^2 # pk is in units of Mpc^3/h^3
   end
   res, err = hquadrature(dσ2dRdk, 0, 20 / Rh)
   return res
end
