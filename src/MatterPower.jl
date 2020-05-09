module MatterPower
using OrdinaryDiffEq
using HCubature
using Roots
export t_nowiggle
export setup_growth
export sigma2, dsigma2dRh
export sigma2gaus, dsigma2gausdRh, d2sigma2gausdRh2
export setup_halofit, halofit
include("eisensteinhu.jl")
include("growth.jl")
include("sigma2.jl")
include("halofit.jl")
end # module
