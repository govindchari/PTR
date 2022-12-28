module PTR
using LinearAlgebra
using ForwardDiff
using Convex, ECOS
using Printf

include("structs.jl")
include("discretize.jl")
include("subproblem.jl")
include("solver.jl")

export ptr, PARAMS, RK4, FOH_discretize!, solveSubproblem!, initialize!, solveTraj!
end