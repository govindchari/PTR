module PTR
using LinearAlgebra
using ForwardDiff
using PyPlot
using Plots
using Convex, ECOS
using Printf

include("structs.jl")
include("discretize.jl")
include("subproblem.jl")
include("solver.jl")
include("postprocess.jl")

export ptr, PARAMS, solveTraj!, plotall, animateTrajectory
end