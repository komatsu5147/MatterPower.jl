module MatterPower
using OrdinaryDiffEqTsit5
using DoubleExponentialFormulas
using Roots
export t_nowiggle
export setup_growth
export sigma2, dsigma2dR
export sigma2gaus, dsigma2gausdR, d2sigma2gausdR2
export setup_halofit, halofit
include("eisensteinhu.jl")
include("growth.jl")
include("sigma2.jl")
include("halofit.jl")
end # module
