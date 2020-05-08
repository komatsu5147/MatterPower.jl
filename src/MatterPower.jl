module MatterPower
using OrdinaryDiffEq
using HCubature
export t_nowiggle, setup_growth
export sigma2, dsigma2dRh
include("eisensteinhu.jl")
include("growth.jl")
include("sigma2.jl")
end # module
