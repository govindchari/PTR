using LinearAlgebra
using Plots
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    # Pendulum Dynamics and Jacobians
    w = 2 * pi
    function f(x::Array{Float64,1}, u::Array{Float64,1})
        return [x[2]; -w^2 * sin(x[1]) + u[1]]
    end
    function dfx(x::Array{Float64,1}, u)
        return [0 1; -w^2*cos(x[1]) 0]
    end
    function dfu(x::Array{Float64,1}, u)
        return [0; 1]
    end

    nx = 2
    nu = 1
    K = 11
    Nsub = 10
    p = PTR.ptr(nx, nu, K, Nsub, f, dfx, dfu)
    PTR.FOH_discretize(p)
    PTR.solveSubproblem(p)

end