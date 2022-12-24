module PTR
using LinearAlgebra
using ForwardDiff
using Convex, ECOS
using Printf

include("structs.jl")
include("discretize.jl")
include("subproblem.jl")

export ptr, RK4, FOH_discretize, solveSubproblem
end